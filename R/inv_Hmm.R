inv_Hmm <- function(A11, A21, A22, G){

  # Assumptions:
  # A11 and A22 are symmetric matrices
  # A21 = t(A12)
  # G is a non-null vector that contains the diagonal elements of the random effects' covariance matrix

  M1 <- try(chol2inv(chol(A11)))
  if (class(M1) == "try-error") {
    M1 <- MASS::ginv(A11)
  }
  M2 <- tcrossprod(A21, M1)
  M3 <- try(chol2inv(chol(diag(1/G) + A22 - tcrossprod(M2, A21))))
  if (class(M3) == "try-error") {
    M3 <- MASS::ginv(diag(1/G) + A22 - tcrossprod(M2, A21))
  }

  S22 <- (1/G)*M3
  # or S22 <- (1/G)*chol2inv(chol(diag(1/G) + A22 - A21 %*% M1 %*% A12))
  S21 <- - S22 %*% M2
  # or S21 <- - S22 %*% A21 %*% M1
  S12 <- - M1 %*% crossprod(A21, M3)
  # or S12 <- - M1 %*% t(A21*G) %*% S22
  S11 <- M1 - S12 %*% M2
  # or S11 <- M1 %*% (diag(ncol(A11)) + t(A21*G) %*% S22 %*% A21 %*% M1)
  S <- sclm::bblock2(S11, S12, S21, S22)

  output <- list(S11 = S11, S12 = S12, S21 = S21, S22 = S22, S = S)
  return(output)
}
