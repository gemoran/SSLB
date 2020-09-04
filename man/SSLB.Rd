\name{SSLB}
\alias{SSLB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spike-and-Slab Lasso Biclustering
}
\description{
\code{SSLB} uses a fast, deterministic EM algorithm to find biclusters of the data matrix,
\code{Y}, via sparse factor analysis where both the factors and the loadings are sparse. 
}
\usage{
SSLB(Y, K_init, 
lambda1 = 1, 
lambda0s = c(1, 5, 10, 50, 100, 500, 1000, 10000, 100000, 1000000, 10000000), 
lambda1_tilde = 1, 
lambda0_tildes = c(1, rep(5, length(lambda0s) - 1)), 
IBP = 1, 
a,  
b = 1,  
a_tilde,  
b_tilde = 1,  
alpha, 
d = 0, 
EPSILON = 0.01, 
MAX_ITER = 500)
}
\arguments{
  \item{Y}{Matrix of observations by features \code{N x G}.
}
  \item{K_init}{Inital number of biclusters.}
  \item{lambda1}{Slab variance parameter for loadings, \code{B}. Needs to be less than 
  or equal to \code{lambda0s}. Default is \code{lambda0 = 1}.}
  \item{lambda0s}{Vector \code{L x 1} of spike parameters for loadings, \code{B}. 
  Default is \code{lambda0s = c(1, 5, 10, 50, 100, 500, 1000, 10000, 100000, 1000000, 
  10000000)}. Note that SSLB relies on a sequence of values for \code{lambda0s}. We do
  NOT recommend using a single value for \code{lambda0s}.}
  \item{lambda1_tilde}{Slab variance parameter for factors, \code{X}. 
  Needs to be less than or equal to \code{lambda0_tildes}. Default is 
  \code{lambda0_tilde = 1}.}
  \item{lambda0_tildes}{Vector \code{L x 1} of spike parameters for factors, \code{X}.
  Default is \code{lambda0_tildes = c(1, rep(5, length(lambda0s) - 1))}.}
  \item{IBP}{Whether an Indian Buffet Process (IBP) prior is used. Either \code{IBP = 1}
  (with IBP prior) or \code{IBP = 0} (with Beta-binomial prior).}
  \item{a}{Hyperparameter for the beta prior \code{B(a, b)} for the loadings, \code{B}.
  Default is \code{a = 1/K_init}.}
  \item{b}{Hyperparameter for the beta prior \code{B(a, b)} for the loadings, \code{B}.}
  \item{a_tilde}{Hyperparameter for the beta prior \code{B(a_tilde, b_tilde)} for the 
  factors, \code{X}. Default is \code{a_tilde = 1/K_init}.}
  \item{b_tilde}{Hyperparameter for the beta prior \code{B(a_tilde, b_tilde)} for the 
  factors, \code{X}.}
 \item{alpha}{Hyperparameter for the IBP prior. Default is \code{alpha = 1/N}.}
 \item{d}{Hyperparameter for the Pitman-Yor IBP prior. Default is \code{d = 0} 
 (reduces to usual IBP).}
 \item{EPSILON}{Convergence criterion: converged when max difference in loadings matrix 
 is less than \code{EPSILON * B_old} (default \code{EPSILON = 0.01}).}
 \item{MAX_ITER}{Maximum number of iterations for each \code{lambda0}.  Default is 500.}
}
\details{
SSLB finds factorization \code{Y = X * t(B)} where both \code{X} and \code{B} are 
column-sparse. This is achieved via a fast EM algorithm. 
}
\value{
A list containing:
\item{X}{The estimated factor matrix \code{X} at the final \code{lambda0s} value 
provided by the user.}
\item{Gamma_tilde}{The estimated support of the factor matrix \code{X} at the final 
\code{lambda0s} value provided by the user.}
\item{B}{The estimated loadings matrix \code{B} at the final \code{lambda0s} value 
provided by the user.}
\item{K}{The estimated number of biclusters at the final \code{lambda0s} value 
provided by the user.}
\item{ML}{The estimated Cholesky decomposition of the variance of \code{X} at the 
final \code{lambda0s} value provided by the user.}
\item{B_init}{Initial \code{B}, generated as 
\code{matrix(rexp(G * K_init, rate = 1), nrow = G, ncol = K_init)}.}
\item{path}{A list of lists containing outputs of the parameters at each value 
of \code{lambda0s}. Parameters are: \code{X}, \code{B}, \code{K}, \code{Tau},
\code{Gamma_tilde}, \code{ML}, \code{theta}, \code{theta_tildes}, \code{nus}, 
\code{sigmas}, \code{lambda0s}, \code{lambda1}, \code{lambda0_tildes}, 
\code{lambda1_tilde}, \code{iter}.}

}
\references{
Moran, G., Rockova, V. and George, E.I. (2019) "Spike-and-Slab Lasso Biclustering" 
Annals of Applied Statistics (Accepted)
}
\author{
Gemma Moran <gm2918@columbia.edu>
}

