phi_inv <- function(z) {
  p <- c(exp(z), 1)
  p / sum(p)
}

lnm_simulator <- function(X, beta, sigma = 1, depth = 3000) {
  N <- nrow(X)
  K <- nrow(beta)
  Z <- rnorm(K, X %*% t(beta), sigma)
  y <- matrix(nrow = N, ncol = K + 1)
  
  for (i in 1:N) {
    z <- rnorm(K, X[i, ] %*% t(beta), sigma)
    y[i, ] <- rmultinom(1, depth, phi_inv(z))
  }
  y
}