h_transform <- function(X, A){
  d <- dim(A)
  M <- matrix(A, nrow = d[1])
  XM <- as.matrix(X %*% M)
  array(XM, c(nrow(XM), d[-1]))
}
