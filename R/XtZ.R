XtZ <- function(X, Z, w = NULL){
  d <- length(Z)
  list.mat <- vector("list", d)
  for (i in d:1) {
    list.mat[[i]] <- sclm::A2_form(X, Z[[i]], w)
    cat("size of block", i, "in XtZ function:", round(object.size(list.mat[[i]])/1048600, 3),"Mb", "\n")
  }
  res <- do.call("cbind", list.mat)
  return(res)
}
