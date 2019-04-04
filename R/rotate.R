rotate <- function(A) {
  d <- 1:length(dim(A))
  d1 <- c(d[-1], d[1])
  aperm(A, d1)
}
