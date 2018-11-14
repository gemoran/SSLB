#ifndef SSLB_H
#define SSLB_H

#include <RcppArmadillo.h>
#include <math.h> 
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

void E_step_X(mat &X, mat &V_traces, mat &ML, mat &B, vec &sigmas, mat &Tau, mat &Y);
mat exp_density(mat &Tau, double lambda);
vec M_step_thetas(mat &Gs, double a, double b);
mat E_step_Gamma_tilde(mat &Tau, vec &theta_tildes, double lambda1_tilde, double lambda0_tilde);
mat M_step_Tau(mat &X, mat &Gamma_tilde, mat &V_traces, double lambda1_tilde, double lambda0_tilde);
void rescale_X_B(mat &X, mat &B);
int check_convergence(mat &B, mat &B_old, double EPSILON) ;

mat E_step_Gamma(mat &B, vec &thetas, double lambda1, double lambda0);

mat M_step_B(mat &Y, mat &B_old, mat &X, mat &ML, vec &sigmas, vec &thetas, double lambda1, double lambda0);
vec SSLASSO(vec &y, mat &X, vec &X_sumsq, vec &beta_old, double lambda1, double lambda0, vec &thetas, double sigma2);
double threshold(double theta, double sigma2, double lambda1, double lambda0, double sumsq);
double g(double x, double theta, double sigma2, double lambda1, double lambda0, double sumsq);
double SSL(double z, double beta, double lambda1, double lambda0, double theta, double sumsq, double delta, double sigma2);
double lambdastar(double x, double theta, double lambda1, double lambda0);
double pstar(double x, double theta, double lambda1, double lambda0);
double sumsq(vec y);
vec M_step_thetas_approx(mat &B, double a, double b);
vec M_step_sigmas(mat &Y, mat &B, mat &X, double nu, double xi, double sigma_min);
vec M_step_nus_IBP(vec &Gamma_tilde_sum, mat &PQ, double alpha, double d, int N);
mat E_step_Q(vec &nus);


#endif
