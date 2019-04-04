Ztheta <- function(Z, theta, np){
  d <- length(Z)
  for (i in 1:d) {
    if (i == 1) {
      res <- sclm::Xtheta(Z[[i]], theta[1:(np[1])])
    } else {
      init <- sum(np[1:(i-1)])
      fin  <- np[i]
      if (fin) res <- res + sclm::Xtheta(Z[[i]], theta[(init+1):(init+fin)])
    }
  }
  res
}
