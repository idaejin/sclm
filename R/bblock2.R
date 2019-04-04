bblock2 <- function(A11, A12, A21, A22){
  block <- rbind(cbind(A11, A12), cbind(A21, A22))
  unname(block, force = FALSE)
}
