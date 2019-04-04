sclm_mat <- function(C, gamma, X, Z, z, w) {
  Ctlist <- lapply(C, Matrix::t)

  CGX. <- sclm::A2_form(Ctlist, X, gamma)
  cat("size of CGX. in sclm_mat function:", round(object.size(CGX.)/1048600, 3), "Mb\n")
  cat("dim of CGX.:", dim(CGX.), "\n")

  CGZ. <- sclm::XtZ(Ctlist, Z, gamma)
  cat("size of CGZ. in clmm_mat_array function:", round(object.size(CGZ.)/1048600, 3), "Mb\n")
  cat("dim of CGZ.:", dim(CGZ.), "\n")

  XtX. <- crossprod(sqrt(1/w)*CGX.) #crossprod(CGX., (1/w)*CGX.)
  XtZ. <- crossprod(CGX., (1/w)*CGZ.)
  ZtX. <- t(XtZ.)
  ZtZ. <- crossprod(sqrt(1/w)*CGZ.)#crossprod(CGZ., (1/w)*CGZ.)
  #ZtZ. <- rfunctions::crossprodcpp(CGZ., 1/w)
  Xty. <- crossprod(CGX., z)
  Zty. <- crossprod(CGZ., z)
  yty. <- sum((z^2)*w)
  ZtXtZ <- rbind(XtZ., ZtZ.)
  u <- c(Xty., Zty.)

  output <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
  return(output)
}
