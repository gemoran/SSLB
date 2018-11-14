\name{SSLB}
\alias{SSLB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spike-and-Slab Lasso Biclustering
}
\description{
\code{SSLB} uses a fast, deterministic EM algorithm to find biclusters of the data matrix, \code{Y},
via sparse factor analysis where both the factors and the loadings are sparse. 

}
\usage{
SSLB(Y, K_init, lambda1, lambda0, lambda1_tilde, lambda0_tildes, IBP = 0,
     a, b = 1, a_tilde, b_tilde = 1, alpha, d = 0.5, sigdf = 1, sigquant = 0.5,
     EPSILON = 0.05, MAX_ITER = 500)
}
\arguments{
  \item{Y}{Matrix of observations by features \code{N x G}.
}
  \item{K_init}{Inital number of biclusters.}
  \item{lambda1}{Slab variance parameter for loadings, \code{B}. Needs to be less than or equal to
  \code{lambda0}. Default is \code{lambda0 = 1}.}
  \item{lambda0}{Vector \code{L x 1} of spike parameters for loadings, \code{B}. Either a numeric value for a single run (L=1) or a sequence of increasing values for dynamic posterior exploration.}
  \item{lambda1_tilde}{Slab variance parameter for factors, \code{X}. Needs to be less than or equal to
  \code{lambda0_tilde}. Default is \code{lambda0_tilde = 1}.}
  \item{IBP}{Whether an Indian Buffet Process (IBP) prior is used. Either \code{IBP = 1} (with IBP prior) or \code{IBP = 0} (with Beta-binomial prior).}
  \item{a}{Hyperparameter for the beta prior \code{B(a, b)} for the loadings, \code{B}.}
  \item{b}{Hyperparameter for the beta prior \code{B(a, b)} for the loadings, \code{B}.}
  \item{a_tilde}{Hyperparameter for the beta prior \code{B(a_tilde, b_tilde)} for the factors, \code{X}.}
  \item{b_tilde}{{Hyperparameter for the beta prior \code{B(a_tilde, b_tilde)} for the factors, \code{X}.}}
 \item{alpha}{}
 \item{d}{}
 \item{sigdf}{}
 \item{sigquant}{}
 \item{EPSILON}{}
 \item{MAX_ITER}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
