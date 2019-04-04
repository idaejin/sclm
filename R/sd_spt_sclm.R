sd_spt_sclm <- function(X, Z, C, Ginv, efine, eta){
  gamma <- c(efine*exp(eta))
  mu <- sclm::Xtheta(C, gamma)
  CGX <- sclm::A2_form(l1 = lapply(C, Matrix::t), l2 = X, w = gamma)
  CGZ <- sclm::XtZ(X = lapply(C, Matrix::t), Z = Z, w = gamma)
  XtX <- crossprod(sqrt(1/mu)*CGX)#crossprod(CGX, (1/mu)*CGX)
  XtZ <- crossprod(CGX, (1/mu)*CGZ)
  ZtX <- t(XtZ)
  ZtZ <- crossprod(sqrt(1/mu)*CGZ)#crossprod(CGZ, (1/mu)*CGZ)
  cov.parts <- sclm::inv_bblock2(XtX, XtZ, ZtX, ZtZ + diag(Ginv))
  S11 <- cov.parts$S11
  S12 <- cov.parts$S12
  S22 <- cov.parts$S22

  m <- length(eta)
  Xt <- lapply(X, t)
  Zt <- vector("list", 7)
  for (j in 1:7) {
    Zt[[j]] <- lapply(Z[[j]], t)
  }

  print(ls())
  rm(gamma, mu, CGX, CGZ, XtX, XtZ, ZtX, ZtZ, C, Ginv, efine, cov.parts, Z)
  gc()
  print(ls())

  var.trend <- rep(0, m)
  vec <- rep(0, m)

  for (i in 1:m) {
    print(i)
    vec[i] <- 1
    Xte <- sclm::Xtheta(Xt, vec)
    Zte.aux <- vector("list", 7)
    for (j in 1:7){
      Zte.aux[[j]] <- sclm::Xtheta(Zt[[j]], vec)
    }
    Zte <- do.call("c", Zte.aux)
    pI <- crossprod(Xte, S11 %*% Xte)
    pII <- crossprod(Xte, S12 %*% Zte)
    pIII <- crossprod(Zte, S22 %*% Zte)
    var.trend[i] <- pI + 2*pII + pIII
    vec[i] <- 0
  }

  sd.trend <- sqrt(var.trend)
  sd.trend.mat <- matrix(sd.trend, nrow(X$X2X1), nrow(X$X3))
  sd.exp.trend <- sd.trend*exp(eta)
  sd.exp.trend.mat <- matrix(sd.exp.trend, nrow(X$X2X1), nrow(X$X3))

  output <- list(sd.trend = sd.trend, sd.trend.mat = sd.trend.mat, sd.exp.trend = sd.exp.trend, sd.exp.trend.mat = sd.exp.trend.mat)
  return(output)
}
