mm_bases <- function(x, xl, xr, ndx, bdeg, pord){
  Bb <- sclm::bspline(x, xl, xr, ndx, bdeg)
  knots <- Bb$knots
  B <- Bb$B
  m <- ncol(B) # number of columns of B
  D <- diff(diag(m), differences = pord)
  P.svd <- svd(crossprod(D)) # equivalent to svd(t(D)%*%D)
  Us <- (P.svd$u)[, 1:(m - pord)] # non-null eigenvectors
  Un <- (P.svd$u)[, -(1:(m - pord))] # null eigenvectors
  d <- (P.svd$d)[1:(m - pord)] # vector of non-null eigenvalues
  Z <- B %*% Us # marginal random model matrix
  X <- B %*% Un  # marginal fixed model matrix

  output <- list(X = X, Z = Z, d = d, B = B, knots = knots, m = m, D = D)
  return(output)
}
