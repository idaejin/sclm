rten2 <- function(X1, X2, sparse = TRUE){
  if (sparse) {
    vec1 <- matrix(1, 1, ncol(X1))
    vec2 <- Matrix::Matrix(1, 1, ncol(X2))
    kronecker(X1, vec2)*Matrix::kronecker(vec1, X2)
  } else {
    vec1 <- matrix(1, 1, ncol(X1))
    vec2 <- matrix(1, 1, ncol(X2))
    kronecker(X1, vec2)*kronecker(vec1, X2)
  }
}
