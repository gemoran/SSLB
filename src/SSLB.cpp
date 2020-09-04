#include "SSLB.h"

// [[Rcpp::export(.cSSLB)]]
SEXP cSSLB(
  SEXP Y_SEXP, 
  SEXP B_init,
  SEXP sigmas_init,
  SEXP Tau_init,
  SEXP thetas_init,
  SEXP theta_tildes_init,
  SEXP nus_init,
  SEXP lambda1_SEXP, 
  SEXP lambda0s_SEXP,
  SEXP lambda1_tilde_SEXP,
  SEXP lambda0_tildes_SEXP,
  SEXP a_SEXP,
  SEXP b_SEXP,
  SEXP a_tilde_SEXP,
  SEXP b_tilde_SEXP,
  SEXP alpha_SEXP,
  SEXP d_SEXP,
  SEXP eta_SEXP,
  SEXP xi_SEXP,
  SEXP sigma_min_SEXP,
  SEXP IBP_SEXP,
  SEXP EPSILON_SEXP,
  SEXP MAX_ITER_SEXP) {

  // Convert R objects to Armadillo objects
  mat Y = as<mat>(Y_SEXP);
  int N = Y.n_rows;
  int G = Y.n_cols;

  mat B = as<mat>(B_init);
  int K_init = B.n_cols;
  int K = K_init;

  mat B_old = B;

  mat X = zeros<mat>(N, K);
  mat V_traces = zeros<mat>(N, K);
  mat ML = zeros<mat>(K, K);

  mat Tau = as<mat>(Tau_init);
  vec sigmas = as<vec>(sigmas_init);
  vec thetas = as<vec>(thetas_init);
  vec nus = as<vec>(nus_init);
  vec theta_tildes = as<vec>(theta_tildes_init); 
  vec lambda0s = as<vec>(lambda0s_SEXP);
  vec lambda0_tildes = as<vec>(lambda0_tildes_SEXP);
  
  mat Gamma_tilde(N, K);
  Gamma_tilde.fill(0.5);

  vec Gamma_tilde_sum(K);
  vec one_to_K = linspace<vec>(0, K - 1, K);

  double lambda1 = as<double>(lambda1_SEXP);
  double lambda1_tilde = as<double>(lambda1_tilde_SEXP);
  double lambda0;
  double lambda0_tilde;
  double a = as<double>(a_SEXP);
  double b = as<double>(b_SEXP);
  double a_tilde = as<double>(a_tilde_SEXP);
  double b_tilde = as<double>(b_tilde_SEXP);
  double alpha = as<double>(alpha_SEXP);
  double d = as<double>(d_SEXP);
  double eta = as<double>(eta_SEXP);
  double xi = as<double>(xi_SEXP);
  double sigma_min = as<double>(sigma_min_SEXP);
  int IBP = as<int>(IBP_SEXP);
  if (IBP == 1) {
    theta_tildes(0) = nus(0);
    for (int k = 1; k < K; k++) {
      theta_tildes(k) = theta_tildes(k - 1) * nus(k);
    }
  }
  mat PQ = zeros<mat>(K, K);

  const int L = lambda0s.n_elem;

  double EPSILON = as<double>(EPSILON_SEXP);
  double MAX_ITER = as<double>(MAX_ITER_SEXP);
  int converged = 0;

  vec NITERS = zeros<vec>(L);
  int ITER = 0;
  int update_lambda0_tilde = 1;

  List B_all(L);
  List X_all(L);
  List ML_all(L);
  List Tau_all(L);
  List Gamma_tilde_all(L);
  mat sigmas_all = zeros<mat>(G, L);
  List thetas_all(L);
  List theta_tildes_all(L);
  List nus_all(L);
  vec Gamma_tilde_sum_all(L);
  vec lambda0_tildes_out(L);

  for(int l = 0; l < L; l++) {

    lambda0 = lambda0s(l);
    if (update_lambda0_tilde == 1) {
      lambda0_tilde = lambda0_tildes(l);
    } 

    Rcout << "Lambda............................................" << l + 1 << endl;

    while(NITERS(l) < MAX_ITER) {
      Rcpp::checkUserInterrupt();

      NITERS(l)++;
      ITER = NITERS(l);


      // Update Gamma_tildes
      if (lambda1 != lambda0) {
        Gamma_tilde = E_step_Gamma_tilde(Tau, theta_tildes, lambda1_tilde, lambda0_tilde);
      }
      
      // Re-order in descending order of Gamma_tilde_sum (for IBP)
      for (int k = 0; k < K; k++) {
        Gamma_tilde_sum(k) = 0;
        for (int i = 0; i < N; i++) {
          Gamma_tilde_sum(k) += Gamma_tilde(i, k);
        }
      }

      if (IBP == 1) {

        uvec Gamma_tilde_order = stable_sort_index(Gamma_tilde_sum, "descend");

        if (any(Gamma_tilde_order != one_to_K)) {
          Gamma_tilde = Gamma_tilde.cols(Gamma_tilde_order);
          Gamma_tilde_sum = Gamma_tilde_sum.elem(Gamma_tilde_order);

          Tau = Tau.cols(Gamma_tilde_order);

          B = B.cols(Gamma_tilde_order);
          B_old = B_old.cols(Gamma_tilde_order);

          thetas = thetas.elem(Gamma_tilde_order);
          nus = nus.elem(Gamma_tilde_order);
        }
      }

      // Update X
      E_step_X(X, V_traces, ML, B, sigmas, Tau, Y);


      // Update B and sigmas
      if (ITER == 1 && l == 0) {
        mat B_temp = zeros<mat>(G, K);
        B = M_step_B(Y, B_temp, X, ML, sigmas, thetas, lambda1, lambda0);
        thetas = M_step_thetas_approx(B, a, b);
      } else {
        B = M_step_B(Y, B_old, X, ML, sigmas, thetas, lambda1, lambda0);
        thetas = M_step_thetas_approx(B, a, b);
      }

      // Update sigmas
      sigmas = M_step_sigmas(Y, B, X, eta, xi, sigma_min);


      // Update theta_tildes
      if (IBP == 1 && lambda1 != lambda0) {
        PQ = E_step_Q(nus);
        nus = M_step_nus_IBP(Gamma_tilde_sum, PQ, alpha, d, N);
        theta_tildes(0) = nus(0);
        for (int k = 1; k < K; k++) {
          theta_tildes(k) = theta_tildes(k - 1) * nus(k);
        }
      } else {
        theta_tildes = M_step_thetas(Gamma_tilde, a_tilde, b_tilde);
      }

      // Update Tau
      Tau = M_step_Tau(X, Gamma_tilde, V_traces, lambda1_tilde, lambda0_tilde);

       // remove zeros
      if ((l != 0 && lambda1 != lambda0) && ITER % 100 == 0) {
        vec B_zero = zeros<vec>(K);
        vec X_zero = zeros<vec>(K);
        for (int k = 0; k < K; k++) {
          for (int j = 0; j < G; j++) {
            if (B(j, k) == 0) {
              B_zero(k)++;
            }
          }
          for (int i = 0; i < N; i++) {
            if (Gamma_tilde(i, k) < 0.025) {
              X_zero(k)++;
            }
          }
        }

        uvec keep = find(B_zero < G-1 && X_zero < N-1);

        if (keep.n_elem < K && keep.n_elem > 0) {
          mat X_new = X.cols(keep);
          K = X_new.n_cols;
          mat Gamma_tilde_new = Gamma_tilde.cols(keep);
          mat Tau_new = Tau.cols(keep);
          mat B_new = B.cols(keep);
          mat B_old_new = B_old.cols(keep);

          Gamma_tilde_sum = Gamma_tilde_sum.elem(keep);
          thetas = thetas.elem(keep);
          theta_tildes = theta_tildes.elem(keep);
          nus = nus.elem(keep);
          uvec new_K = linspace<uvec>(0, K - 1, K);
          one_to_K = one_to_K.elem(new_K);

          V_traces.set_size(N, K);
          ML.set_size(K, K);
          PQ.set_size(K, K);
          X.set_size(size(X_new));
          X = X_new;
          Tau.set_size(size(Tau_new));
          Tau = Tau_new;
          B.set_size(size(B_new));
          B = B_new;
          B_old.set_size(size(B_old_new));
          B_old = B_old_new;
          Gamma_tilde.set_size(size(Gamma_tilde_new));
          Gamma_tilde = Gamma_tilde_new;
        }

        if (keep.n_elem == 0) {
          K = 0;
          break;
        }

      }
      
      // remove zero components
      if (K == 0) {
        Rcout << "Number of biclusters is 0" << endl;
        break;
      }

      // Check convergence
      converged = check_convergence(B, B_old, EPSILON);

      if (converged == 1) {
        Rcout << "Iteration: " << ITER  << "......Converged." << endl;
        break;
      }

      B_old = B;

      if (ITER % 100 == 0) {
        Rcout << "Iteration: " << ITER << endl;

      }

      // Rescale X and B 
     if (lambda0 != lambda1) {
      rescale_X_B(X, B);
     }
      
    }

    lambda0_tildes_out(l) = lambda0_tilde;

    Gamma_tilde_sum_all(l) = sum(Gamma_tilde_sum);
    if (l > 1 && update_lambda0_tilde == 1) {
      if (Gamma_tilde_sum_all(l) > Gamma_tilde_sum_all(l - 1)) {
        lambda0_tilde = lambda0_tildes(l - 1);
        update_lambda0_tilde = 0;
      }
    }

    B_all[l] = B;
    X_all[l] = X;
    ML_all[l] = ML;
    Tau_all[l] = Tau;
    Gamma_tilde_all[l] = Gamma_tilde;
    sigmas_all.col(l) = sigmas;
    thetas_all[l] = thetas;
    theta_tildes_all[l] = theta_tildes;
    nus_all[l] = nus;

    if (K == 0) {
      break;
    }

    // if (ITER <= 5) {
    //   break;
    // }

  }

  List out(20);
  out["X"] = X_all;
  out["B"] = B_all;
  out["ML"] = ML_all;
  out["Tau"] = Tau_all;
  out["Gamma_tilde"] = Gamma_tilde_all;
  out["sigmas"] = sigmas_all;
  out["thetas"] = thetas_all;
  out["theta_tildes"] = theta_tildes_all;
  out["nus"] = nus_all;
  out["lambda1"] = lambda1;
  out["lambda0s"] = lambda0s;
  out["lambda1_tilde"] = lambda1_tilde;
  out["lambda0_tildes"] = lambda0_tildes_out;
  out["a"] = a;
  out["b"] = b;
  out["a_tilde"] = a_tilde;
  out["b_tilde"] = b_tilde;
  out["nu"];
  out["xi"];
  out["iter"] = NITERS;

  return out;

}