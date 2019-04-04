A2_form <- function(l1, l2, w = NULL){

  d1 <- length(l1)
  d2 <- length(l2)

  if(!(d1 == d2)) {
    stop("l1 and l2 should have the same dimension")
  }

  n <- rev(sapply(l1, nrow))
  d <- rev(sapply(l1, ncol))
  c <- rev(sapply(l2, ncol))

  if (is.null(w)) {
    W <- array(1, n)
  } else {
    W <- array(w, n)
  }

  tmp <- sclm::rh(Matrix::t(sclm::rten2(l2[[d1]], l1[[d1]])), W)

  for (i in (d1-1):1) {
    tmp <- sclm::rh(Matrix::t(sclm::rten2(l2[[i]], l1[[i]])), tmp)
  }

  dim(tmp)<- as.vector(rbind(d,c))
  Fast1 <- aperm(tmp, as.vector(matrix(1:(d1*2), byrow = TRUE, ncol = 2)))

  Fast <- if(prod(d)) matrix(Fast1, nrow = prod(d)) else Fast1
  return(Fast)
}
