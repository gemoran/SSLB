
SSLB <- function(Y, 
                 K_init,
                 lambda1,
                 lambda0s,
                 lambda1_tilde,
                 lambda0_tildes,
                 IBP = 1,
                 a, 
                 b = 1, 
                 a_tilde, 
                 b_tilde = 1, 
                 alpha,
                 d = 0,
                 EPSILON = 0.01, 
                 MAX_ITER = 500) {
  
  N <- nrow(Y)
  G <- ncol(Y)
  
  if (missing(K_init)) {
    stop("Must provide initial value of K (K_init)")
  }

  if (missing(a)) {
    a <- 1/(K_init)
  }
  if (missing(a_tilde)) {
    a_tilde <- 1/(K_init)
  }

  if (missing(alpha)) {
    alpha <- 1/N
  }
  
  sigs <- apply(Y, 2, sd)
  
  sigquant <- 0.5
  sigdf <- 3
  
  sigest <- quantile(sigs, 0.05)
  qchi <- qchisq(1 - sigquant, sigdf)
  xi <- sigest^2 * qchi / sigdf
  eta <- sigdf
  sigmas_median <- sigest^2
  sigmas_init <- rep(sigmas_median, G)
  sigma_min <- sigest^2 / G
  
  
  B_init <- matrix(rexp(G * K_init, rate = 1), nrow = G, ncol = K_init)
  Tau_init <- matrix(100, nrow = N, ncol = K_init)
  thetas_init <- rep(0.5, K_init)
  nus_init <- sort(rbeta(K_init, 1, 1), decreasing = T)
#  nus_init <- c(0.5, rep(1, K_init-1))
  theta_tildes_init <- rep(0.5, K_init)
  
  nlambda <- length(lambda0s)
  
  res <- .Call("cSSLB", Y, B_init, sigmas_init, Tau_init, thetas_init, theta_tildes_init, 
    nus_init, lambda1, lambda0s, lambda1_tilde, lambda0_tildes, a, b, 
        a_tilde, b_tilde, alpha, d, eta, xi, sigma_min, IBP, EPSILON, MAX_ITER, PACKAGE = "SSLB")
  
  X <- res$X
  X <- X[!sapply(X, is.null)]
  Tau <- res$Tau
  Tau <- Tau[!sapply(Tau, is.null)]
  Gamma_tilde <- res$Gamma_tilde
  Gamma_tilde <- Gamma_tilde[!sapply(Gamma_tilde, is.null)]
  B <- res$B
  B <- B[!sapply(B, is.null)]
  ML <- res$ML
  ML <- ML[!sapply(ML, is.null)]
  
  thetas <- lapply(res$thetas, as.vector)
  thetas <- thetas[!sapply(thetas, is.null)]
  
  theta_tildes <- lapply(res$theta_tildes, as.vector)
  theta_tildes <- theta_tildes[!sapply(theta_tildes, is.null)]
  
  nus <- lapply(res$nus, as.vector)
  nus <- nus[!sapply(nus, is.null)]
  
  sigmas <- res$sigmas
  
  iter <- as.vector(res$iter)
  
  keep_path <- which(iter > 0)
  
  lambda0s <- lambda0s[keep_path]
  lambda0_tildes <- res$lambda0_tildes[keep_path]
  iter <- iter[keep_path]
  sigmas <- sigmas[, keep_path]
  nlambda <- length(lambda0s)
  
  K <- numeric(nlambda)
  
  for (l in 1:nlambda) {
    if (lambda1 != lambda0s[l]) {
      X[[l]][Gamma_tilde[[l]] < 0.5] <- 0
      Tau[[l]][Gamma_tilde[[l]] < 0.5] <- 0
      keep <- which(apply(X[[l]], 2, function(x) sum(x != 0) > 1))
      X[[l]] <- as.matrix(X[[l]][, keep])
      B[[l]] <- as.matrix(B[[l]][, keep])
      ML[[l]] <- as.matrix(ML[[l]][keep, keep])
      Tau[[l]] <- as.matrix(Tau[[l]][, keep])
      Gamma_tilde[[l]] <- as.matrix(Gamma_tilde[[l]][, keep])
      thetas[[l]] <- thetas[[l]][keep]
      theta_tildes[[l]] <- theta_tildes[[l]][keep]
      nus[[l]] <- nus[[l]][keep]
      
      K[l] <- length(keep)
    }

  }
  
  
  
  names_out <- c("X", "B", "Gamma_tilde", "ML", "K", "path", "init_B")
  names_path <- c("X", "B", "K", "Tau", "Gamma_tilde", "ML", "thetas", 
                  "theta_tildes", "sigmas", "lambda0s", "lambda1", 
                  "lambda0_tildes", "lambda1_tilde", "iter")
  
  out <- vector("list", length(names_out))
  names(out) <- names_out
  path <- vector("list", length(names_path))
  names(path) <- names_path
  path$X <- X
  path$B <- B
  path$K <- K
  path$Tau <- Tau
  path$Gamma_tilde <- Gamma_tilde
  path$ML <- ML
  path$thetas <- thetas
  path$theta_tildes <- theta_tildes
  path$nus <- nus
  path$sigmas <- sigmas
  path$lambda0s <- lambda0s
  path$lambda1 <- lambda1
  path$lambda0_tildes <- lambda0_tildes
  path$lambda1_tilde <- lambda1_tilde
  path$iter <- iter
  
  out$path <- path
  out$X <- X[[nlambda]]
  out$Gamma_tilde <- Gamma_tilde[[nlambda]]
  out$B <- B[[nlambda]]
  out$K <- K[nlambda]
  out$ML <- ML[[nlambda]]
  out$init_B <- B_init

  return(out)
}