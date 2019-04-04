bspline <- function(x, xl, xr, ndx, bdeg){
  dx <- (xr - xl)/ndx
  knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by = dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1)*dx^bdeg)
  B <- (-1)^(bdeg + 1)*P %*% t(D)

  output <- list(B = B, knots = knots)
  return(output)
}
