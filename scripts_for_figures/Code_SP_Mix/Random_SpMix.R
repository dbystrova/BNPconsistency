##########################################################################
##Generating observations from a multivariate normal mixture distribution:
##########################################################################

MultVar_RandomNormal <- function(N, eta, mu, sigma, z,seed_ =1000) {
  set.seed(seed_)
  K <- length(mu[1, ])
  r <- length(mu[, 1])
  y <- matrix(0, N, r)
  if (sum(z) == FALSE) {
    z <- sample(1:K, N, replace = TRUE, prob = eta)
  }
  for (k in 1:K) {
    if (sum(z == k)) {
      y[z == k, ] <- mvrnorm(sum(z == k), mu[, k], sigma[, , k])
    }
  }
  return(list(y = y, z = z, eta = eta, mu = mu, sigma = sigma))
}



