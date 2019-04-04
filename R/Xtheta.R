Xtheta <- function(X, theta){
  d <- length(X)
  n <- rev(sapply(X, ncol))
  Theta <- array(theta, n)
  tmp <- sclm::rh(X[[d]], Theta)
  for (i in (d-1):1) {
    tmp <- sclm::rh(X[[i]], tmp)
  }
  as.vector(tmp)
}
