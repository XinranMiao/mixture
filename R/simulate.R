phi_inv <- function(mu) {
  p <- c(exp(mu), 1)
  p / sum(p)
}

lnm_simulator <- function(X, beta, sigma = 1, depth = 3000) {
  N <- nrow(X)
  K <- nrow(beta)
  mu <- rnorm(K, X %*% t(beta), sigma)
  y <- matrix(nrow = N, ncol = K + 1)
  
  for (i in 1:N) {
    mu <- rnorm(K, X[i, ] %*% t(beta), sigma)
    y[i, ] <- rmultinom(1, depth, phi_inv(mu))
  }
  y
}