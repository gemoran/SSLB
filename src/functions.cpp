#include "SSLB.h"

// check convergence

int check_convergence(mat &B, mat &B_old, double EPSILON) {
  int G = B.n_rows;
  int K = B.n_cols;
  double diff = 0;
  int converged = 1;
  for (int j = 0; j < G; j++) {
    for (int k = 0; k < K; k++) {
      diff = fabs((B(j, k) - B_old(j, k))/B_old(j, k));
      if (diff > EPSILON) {
        converged = 0;
        break;
      }
    }
    if (converged == 0) {
      break;
    }
  }
  return converged;
}

// E-Step for factor matrix X
void E_step_X(mat &X, mat &V_traces, mat &ML, mat &B, vec &sigmas, mat &Tau, mat &Y) { 
  
  int K = Tau.n_cols;
  mat SinvB = B.each_col() / sigmas;
  mat tBSinvB = trans(B) * SinvB;
  mat M = zeros<mat>(K, K);
  mat V = zeros<mat>(K, K);
  mat Vinv = zeros<mat>(K, K);
  vec X_mean = zeros<vec>(K);
  vec e(K);
  e.fill(1e-5);
  mat keep_sympd = diagmat(e);

  for(int i = 0; i < Tau.n_rows; i++) {
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }

    Vinv = tBSinvB + diagmat(1/Tau.row(i));
    V = inv_sympd(Vinv + keep_sympd);
    V_traces.row(i) = trans(diagvec(V));
    X_mean = V * trans(B) * (1/sigmas % trans(Y.row(i)));
    M += V;
    X.row(i) = trans(X_mean);
  }

  ML = chol(M, "lower");
}

// E-Step for loadings indicator matrix Gamma
mat E_step_Gamma(mat &B, vec &thetas, double lambda1, double lambda0) {
  mat mult_mat = zeros<mat>(B.n_rows, B.n_cols);
  mult_mat.each_row() = trans((1 - thetas)/(thetas + 1e-8));
  mult_mat = (lambda0/lambda1) * mult_mat;
  mat denom = 1 + mult_mat % exp(-abs(B) *  (lambda0 - lambda1));
  mat Gs = 1 / denom;
  return Gs;
}

// E-Step for factor indicator matrix Gamma_tilde
mat E_step_Gamma_tilde(mat &Tau, vec &theta_tildes, double lambda1_tilde, double lambda0_tilde) {
  int N = Tau.n_rows;
  int K = Tau.n_cols;
  mat mult_mat = zeros<mat>(N, K);
  mult_mat.each_row() = trans((1 - theta_tildes)/theta_tildes);
  mult_mat = (pow(lambda0_tilde, 2)/pow(lambda1_tilde, 2)) * mult_mat;
  mat denom = 1 + mult_mat % exp(-Tau * (pow(lambda0_tilde, 2) - pow(lambda1_tilde, 2)) / 2);
  mat Gamma_tilde = 1/denom;
  return Gamma_tilde;
}

// M-step for error variances sigmas
vec M_step_sigmas(mat &Y, mat &B, mat &X, double eta, double xi, double sigma_min) {
  int G = Y.n_cols;
  int N = Y.n_rows;
  vec sigmas = zeros<vec>(G);
  double resid = 0;
  for (int j = 0; j < G; j++) {
    resid = sum(pow(Y.col(j) - X * trans(B.row(j)), 2));
    sigmas(j) = (resid + eta * xi) / (N + eta + 2);
    if (sigmas(j) < sigma_min) {
      sigmas(j) = eta * xi / (eta + 2);
    }
  }
  return sigmas;
}

// Functions for SSLASSO
double sumsq(vec y) {
  int n = y.n_elem;
  double sumsq = 0;
  for(int i = 0; i < n; i++) {
    sumsq += pow(y(i), 2);
  }
  return sumsq;
}

double pstar(double x, double theta, double lambda1, double lambda0){
  double value;
  if (lambda1 == lambda0){
    return 1;
  } else{
    value = (1 - theta)/theta * lambda0 / lambda1 * exp(-fabs(x) * (lambda0 - lambda1));
    value += 1;
    value = 1/value;
    return value;
  }
}

double lambdastar(double x, double theta, double lambda1, double lambda0){
  double aux;
  if (lambda1 == lambda0){
    return lambda1;
  } else{
    aux = pstar(x, theta, lambda1, lambda0);
    return aux * lambda1 + (1 - aux) * lambda0;
  }
}

double SSL(double z, double beta, double lambda1, double lambda0, double theta, double sumsq, double delta, double sigma2) {
  double s=0;
  double lambda;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= delta) {
    return(0);
  } else { 
    lambda=lambdastar(beta, theta, lambda1, lambda0);
    double temp;
    temp = fabs(z) - sigma2 * lambda;
    if (temp > 0) {
      return(temp * s / sumsq);
    } else {
      return(0);  
    }
  }
}

double g(double x, double theta, double sigma2, double lambda1, double lambda0, double sumsq){
  double value=lambdastar(x,theta,lambda1,lambda0); 
  return pow((value-lambda1),2)+2*sumsq/sigma2*log(pstar(x,theta,lambda1,lambda0));
}


double threshold(double theta, double sigma2, double lambda1, double lambda0, double sumsq){ 
  if (lambda0 == lambda1){
    return sigma2 * lambda1;
  } else {
    if (g(0, theta, sigma2, lambda1 ,lambda0, sumsq) > 0){
      return sqrt(2 * sumsq * sigma2 * log(1/pstar(0, theta, lambda1, lambda0))) + sigma2 * lambda1;     
    }  else {    
      return sigma2 * lambdastar(0, theta, lambda1, lambda0);     
    }
  }
}

// SSLASSO
vec SSLASSO(vec &y, mat &X, vec &X_sumsq, vec &beta_old, double lambda1, double lambda0, vec &thetas, double sigma2) {
  int K = X.n_cols;
  vec beta = beta_old;
  double eps = 0.05;
  int max_iter = 100;
  int iter = 0;
  vec r = y - X * beta;
  vec z = zeros<vec>(K);
  double delta;

  while (iter < max_iter) {
    iter++;

    for(int k = 0; k < K; k++) {
      delta = threshold(thetas(k), sigma2, lambda1, lambda0, X_sumsq(k));
      z(k) = dot(X.col(k), r) + X_sumsq(k) * beta_old(k);
      beta(k) = SSL(z(k), beta_old(k), lambda1, lambda0, thetas(k), X_sumsq(k), delta, sigma2);
      double shift = beta(k) - beta_old(k);
      if (shift != 0) {
        r -= shift * X.col(k);
      }
    } 
    if (norm(beta - beta_old, 2) < eps) {
      break;
    }
    beta_old = beta;    
  }
  return beta;
}


mat M_step_B(mat &Y, mat &B_old, mat &X, mat &ML, vec &sigmas, vec &thetas, double lambda1, double lambda0) {
  int K = X.n_cols;
  int G = Y.n_cols;
  mat Y_star = join_cols(Y, zeros<mat>(K, G));
  mat X_star = join_cols(X, ML);
  mat B = zeros<mat>(G, K);
  vec beta_old;
  vec X_sumsq(K);
  for(int k = 0; k < K; k++) {
    X_sumsq(k) = sumsq(X_star.col(k));
  }
  vec Yj;
  vec out = zeros<vec>(K);
  for(int j = 0; j < G; j++) {
    if (j % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    Yj = Y_star.col(j);
    beta_old = trans(B_old.row(j));
    out = SSLASSO(Yj, X_star, X_sumsq, beta_old, lambda1, lambda0, thetas, sigmas(j));
    B.row(j) = trans(out);
  }
  return B;
}

vec M_step_thetas_approx(mat &B, double a, double b) {
  int K = B.n_cols;
  int G = B.n_rows;
  vec thetas = zeros<vec>(K);
  vec sum_B = zeros<vec>(K);
  for(int k = 0; k < K; k++){
    for(int j = 0; j < G; j++){
      if(B(j, k) != 0) sum_B(k)++;
    }
  }
  thetas = (a + sum_B)/(a + b + G);
  for (int k = 0; k < K; k++) {
    if (thetas(k) <= 0) {
      thetas(k) = 1e-10;
    } 
    if (thetas(k) > 1) {
      thetas(k) = 0.99999;
    }
  }
  return thetas;
}

vec M_step_thetas(mat &Gs, double a, double b) {
  int K = Gs.n_cols;
  vec thetas = zeros<vec>(K);
  int N = Gs.n_rows;
  double G_sum = 0;
  for (int k = 0; k < K; k++) {
    for (int i = 0; i < N; i++) {
      G_sum += Gs(i, k);
    }
    thetas(k) = (a + G_sum - 1) / (a + b + N - 2);
    if (thetas(k) < 0) {
      thetas(k) = 1e-10;
    } 
    if (thetas(k) > 1) {
      thetas(k) = 0.99999;
    }
    G_sum = 0;
  }
  return thetas;
}

mat M_step_Tau(mat &X, mat &Gamma_tilde, mat &V_traces, double lambda1_tilde, double lambda0_tilde) {
  int n = X.n_rows;
  int k = X.n_cols;
  mat lambda0s = zeros<mat>(n, k);
  lambda0s.fill(pow(lambda0_tilde, 2));
  mat lambda1s = zeros<mat>(n, k);
  lambda1s.fill(pow(lambda1_tilde, 2));
  mat Lambda_star = lambda1s % Gamma_tilde + (1 - Gamma_tilde) % lambda0s; 
  mat Tau = (-1 + sqrt(1 + 4 * Lambda_star % (pow(X, 2) + V_traces)))/(2 * Lambda_star);
  Tau.elem( find_nonfinite(Tau) ).fill(1e-8);
  return Tau;
}


// Rescale X and B 
void rescale_X_B(mat &X, mat &B) {
  int K = X.n_cols;
  int N = X.n_rows;
  int G = B.n_rows;
  double X_norm = 0; 
  double B_norm = 0;
  rowvec d = ones<rowvec>(K);

  for(int k = 0; k < K; k++) {
    for (int i = 0; i < N; i++) {
      X_norm += fabs(X(i, k));
    }
    for (int j = 0; j < G; j++) {
      B_norm += fabs(B(j, k));
    }
    if ((X_norm > 0) && (B_norm > 0)) {
      d(k) = pow(X_norm / B_norm, 0.5);
    }
    X_norm = 0;
    B_norm = 0;
  }
  X.each_row() /= d;
  B.each_row() %= d;

}


mat E_step_Q(vec &nus) {
  int K = nus.n_elem;
  mat PQ = zeros<mat>(K, K);
  PQ(0, 0) = 1;
  double mult = 1;
  for (int k = 1; k < K; k++) {
    PQ(k, 0) = 1 - nus(0);
    mult = 1;
    for (int l = 1; l <= k; l++) {
      mult *= nus(l - 1);
      PQ(k, l) = (1 - nus(l)) * mult;
    }
  }
  colvec PQ_sum = sum(PQ, 1);
  PQ_sum(find(PQ_sum == 0)).ones();
  PQ.each_col() /= PQ_sum;

  return PQ;
}




vec M_step_nus_IBP(vec &Gamma_tilde_sum, mat &PQ, double alpha, double d, int N) {
  int K = Gamma_tilde_sum.n_elem;
  vec nus = zeros<vec>(K);
  double a = 0;
  double b = 0;
  double PQ_sum = 0;

  for (int k = 0; k < K; k++) {
    for (int m = k; m < K; m++) {
      a += Gamma_tilde_sum(m);
      b += (N - Gamma_tilde_sum(m)) * PQ(m, k);
    }
    for (int m = k + 1; m < K; m++) {
      for (int i = k + 1; i <= m; i++) {
        PQ_sum += PQ(m, i);
      }
      a += (N - Gamma_tilde_sum(m)) * PQ_sum;
      PQ_sum = 0;
    }
    a += alpha + k * d - 1;
    b += - d;
    nus(k) = a/(a + b);
    if (nus(k) > 1) {
      nus(k) = (a + 1)/(a + b + 2);
    }
    a = 0;
    b = 0;
    if (nus(k) < 0) {
      nus(k) = 1e-8; 
    }
  }
  uvec replace_nus = find_nonfinite(nus);
  for (int i = 0; i < replace_nus.n_elem; i++) {
    nus(replace_nus(i)) = 1/(replace_nus(i) + 1);
  }

  return nus;
}





// IBP using NLOPT

// static int fcount = 0, ccount = 0;

// typedef struct {
//   std::vector<double> coeffs;
//   int N;
//   double alpha;
// } my_func_data;


// double myfunc(unsigned n, const double* x, double* grad, void* f_data) {
//   fcount++;
//   my_func_data *d = (my_func_data *) f_data;
//   std::vector<double> coeffs = d->coeffs;
//   int N = d->N;
//   double alpha = d->alpha;

//   if (grad) {
//     for (int i = 0; i < n; i++) {
//       grad[i] = -(coeffs[i] - 1)/x[i] + (N - coeffs[i])/(1 - x[i]);
//     }
//     grad[0] -= (alpha) / x[n-1];
//   }
  
//   double fun = 0;
//   for (int i = 0; i < n; i++) {
//     fun += (coeffs[i] - 1) * log(x[i]) + (N - coeffs[i]) * log(1 - x[i]);
//   }
//   fun += (alpha) * log(x[n-1]);
  

//   return -fun;
// }


// void multi_constraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
//       ccount++;

//   // n is the length of x, m is the length of result (that is, m is number of constraints, Ci)

//   // The n dimension of grad is stored contiguously, so that \partci/\partxj is stored in grad[i*n + j]
//   // Here you see take dCi/dx0...dxn and store it one by one, then repeat. grad is just an one dimensional array

//   if (grad) {
//     for (int i = 0; i < m; i++) {
//       grad[i * n + i] = -1;
//       grad[i * n + (i + 1)] = 1;
//     }
//   }
  
//   for (int i = 0; i < m; i++) {
//     result[i] = x[i + 1] - x[i];
//   }

// }

// vec M_step_theta_tildes_IBP(mat Gamma_tilde, vec theta_tildes, double alpha) {

//   int N = Gamma_tilde.n_rows;
//   int K = Gamma_tilde.n_cols;

//   std::vector<double> coeffs(K);

//   for(int k = 0; k < K; k++) {
//     coeffs[k] = 0;
//     for(int i = 0; i < N; i++) {
//       coeffs[k] += Gamma_tilde(i, k);
//     }
//   }

//   nlopt_opt opt;

//   opt = nlopt_create(NLOPT_LN_COBYLA, K); 
  
//   double lb[K];
//   double ub[K];
//   double tol[K];
  
//   for (int k = 0; k < K; k++) {
//     lb[k] = 0;
//     ub[k] = 1;
//     tol[k] = 1e-6;
//   }

//   nlopt_set_lower_bounds(opt, lb);
//   nlopt_set_upper_bounds(opt, ub);
  
//   my_func_data f_data;
//   f_data.coeffs = coeffs;
//   f_data.N = N;
//   f_data.alpha = alpha;
  
//   nlopt_set_min_objective(opt, myfunc, &f_data);
//   nlopt_add_inequality_mconstraint(opt, K - 1, multi_constraint, NULL, tol);

//   nlopt_set_ftol_rel(opt, 1e-6);
//   nlopt_set_maxeval(opt, 500);

//   std::vector<double> x(K);
//   for (int i = 0; i < K; i++) {
//     x[i] = theta_tildes(i);
//   }

//   double minf;
  
//   if (nlopt_optimize(opt, &(x[0]), &minf) < 0) {
//    for (int k = 0; k < K; k++) {
//     theta_tildes(k) = 0.5; 
//   }
// } else {
//   theta_tildes = arma::conv_to<arma::vec>::from(x);

// } 

// nlopt_destroy(opt);

// fcount = 0;
// ccount = 0;

// return theta_tildes;


// }
