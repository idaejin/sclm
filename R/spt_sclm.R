spt_sclm <- function(y, x1, x2, x3, efine = NULL, Cs, Ct, ndx = c(10, 10, 10), bdeg = c(3, 3, 3), pord = c(2, 2, 2), thr = c(1e-06, 1e-06), maxit = c(100, 300), trace = TRUE, bold = NULL){

  ### TO DO list:
  # - Check inputs
  # - Put warnings if necessary.

  # Set initial time
  start.all <- proc.time()[3]

  # Create x1lim, x2lim and x3lim variables
  x1lim <- c(min(x1) - 0.01, max(x1) + 0.01)
  x2lim <- c(min(x2) - 0.01, max(x2) + 0.01)
  x3lim <- c(min(x3) - 0.01, max(x3) + 0.01)

  # Check for efine variable
  dimfine <- length(x1)*length(x3)
  if (is.null(efine)) {
    efine <- rep(1, dimfine)
  } else if (length(efine) == 1) {
    efine <- rep(efine, dimfine)
  } else {
    efine <- efine
  }

  # Setup for mixed model representation
  MM1 <- sclm::mm_bases(x1, x1lim[1], x1lim[2], ndx[1], bdeg[1], pord[1])
  MM2 <- sclm::mm_bases(x2, x2lim[1], x2lim[2], ndx[2], bdeg[2], pord[2])
  MM3 <- sclm::mm_bases(x3, x3lim[1], x3lim[2], ndx[3], bdeg[3], pord[3])

  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; c1 <- MM1$m
  X2 <- MM2$X; Z2 <- MM2$Z;	d2 <- MM2$d; c2 <- MM2$m
  X3 <- MM3$X; Z3 <- MM3$Z;	d3 <- MM3$d; c3 <- MM3$m

  # Create elements for the inverse of the covariance matrix G
  g1u <- rep(1, pord[3])%x%rep(1, pord[2])%x%d1
  g2u <- rep(1, pord[3])%x%d2%x%rep(1, pord[1])
  g3u <- d3%x%rep(1, pord[2])%x%rep(1, pord[1])

  g11b <- rep(1, pord[3])%x%rep(1, c2-pord[2])%x%d1
  g12b <- rep(1, c3-pord[3])%x%rep(1, pord[2])%x%d1

  g21b <- rep(1, pord[3])%x%d2%x%rep(1, c1-pord[1])
  g22b <- rep(1, c3-pord[3])%x%d2%x%rep(1, pord[1])

  g31b <- d3%x%rep(1, pord[2])%x%rep(1, c1-pord[1])
  g32b <- d3%x%rep(1, c2-pord[2])%x%rep(1, pord[1])

  g1t <- rep(1, c3-pord[3])%x%rep(1, c2-pord[2])%x%d1
  g2t <- rep(1, c3-pord[3])%x%d2%x%rep(1, c1-pord[1])
  g3t <- d3%x%rep(1, c2-pord[2])%x%rep(1, c1-pord[1])

  # Number of parameters in each part
  np <- c(prod(pord),
          (c3-pord[3])*prod(pord[1:2]),
          (c2-pord[2])*prod(pord[c(1,3)]),
          (c1-pord[1])*prod(pord[2:3]),
          (c3-pord[3])*(c2-pord[2])*pord[1],
          (c3-pord[3])*(c1-pord[1])*pord[2],
          (c2-pord[2])*(c1-pord[1])*pord[3],
          (c3-pord[3])*(c2-pord[2])*(c1-pord[1]))

  # Create capital lambda matrices
  G1inv.n <- c(rep(0, sum(np[2:3])), g1u, rep(0, np[5]), g12b, g11b, g1t)
  G2inv.n <- c(rep(0, np[2]), g2u, rep(0, np[4]), g22b, rep(0, np[6]), g21b, g2t)
  G3inv.n <- c(g3u, rep(0, sum(np[3:4])), g32b, g31b, rep(0, np[7]), g3t)

  # 0	0
  # 0 I
  #D <- diag(c(rep(0, np[1]), rep(1, sum(np[-1]))))

  # Set initial variance components
  la <- c(1, 1, 1)

  # Set initial fixed and random effects
  if (is.null(bold)) {
    bold <- rep(0, sum(np))
  }

  # Create lists of spatio-temporal mixed model matrices
  X2X1 <- sclm::rten2(X2, X1, sparse = FALSE)
  Z2X1 <- sclm::rten2(Z2, X1, sparse = FALSE)
  X2Z1 <- sclm::rten2(X2, Z1, sparse = FALSE)
  Z2Z1 <- sclm::rten2(Z2, Z1, sparse = FALSE)

  X <- list(X3 = X3, X2X1 = X2X1)
  Z <- list(Z1 = list(Z3 = Z3, X2X1 = X2X1),
            Z2 = list(X3 = X3, Z2X1 = Z2X1),
            Z3 = list(X3 = X3, X2Z1 = X2Z1),
            Z4 = list(Z3 = Z3, Z2X1 = Z2X1),
            Z5 = list(Z3 = Z3, X2Z1 = X2Z1),
            Z6 = list(X3 = X3, Z2Z1 = Z2Z1),
            Z7 = list(Z3 = Z3, Z2Z1 = Z2Z1))

  # Build list of marginal composition matrices
  # *Note: Usually, Cs is sparse and its dimension is bigger dimension than the dimension of Ct
  C <- list(Ct = Ct, Cs = Matrix::Matrix(Cs))

  # Compute the trend at fine scale
  eta <- sclm::Xtheta(X, bold[1:np[1]]) + sclm::Ztheta(Z, bold[-(1:np[1])], np[-1])

  gamma <- c(efine*exp(eta)) # equivalent to c(exp(eta+log(efine)))
  mu <- sclm::Xtheta(C, gamma)

  # Loop for eta ...
  for (i in 1:(maxit[1])) {

    if (trace) { start <- proc.time()[3] }

    geta <- c(gamma*eta)
    z <- (Xtheta(X = C, theta = geta) +  y - mu)/mu
    mat <- sclm::sclm_mat(C, gamma, X, Z, z, mu)

    # Loop for variance components
    for (it in 1:(maxit[2])) {

      # Create (the diagonal of the) block diagonal matrix
      Ginv <- c((1/la[3])*g3u,
                (1/la[2])*g2u,
                (1/la[1])*g1u,
                (1/la[2])*g22b + (1/la[3])*g32b,
                (1/la[1])*g12b + (1/la[3])*g31b,
                (1/la[1])*g11b + (1/la[2])*g21b,
                (1/la[1])*g1t + (1/la[2])*g2t + (1/la[3])*g3t)
      G <- 1/Ginv

      # X'X X'ZG
      # Z'X Z'ZG
      #V <- bblock2(mat$XtX., t(mat$ZtX.*G), mat$ZtX., t(mat$ZtZ.*G))
      # V + D
      #V <- V + D

      # Hinv <- try(solve(V))
      # if (class(Hinv) == "try-error") {
      #   Hinv <- MASS::ginv(V)
      # }
      #print("Hinv was calculated")

      Hinv <- sclm::inv_Hmm(mat$XtX., mat$ZtX., mat$ZtZ., G)$S
      # Compute new fixed and random effects
      b <- Hinv %*% mat$u

      b.fixed <- b[1:np[1]]
      b.random <- G*b[-(1:np[1])]

      # Compute effective dimensions and variance components
      # Only the diagonal of ZtNZ
      dZtNZ <- apply((t(Hinv[-(1:np[1]),])*mat$ZtXtZ), 2, sum)

      # Tau 1
      G1inv.d <- (1/la[1])*G1inv.n
      ed1 <- sum(dZtNZ*(G1inv.d*G^2))
      ed1 <- ifelse(ed1 == 0, 1e-50,ed1)
      tau1 <- sum(b.random^2*G1inv.n)/ed1
      tau1 <- ifelse(tau1 == 0, 1e-50,tau1)

      # Tau 2
      G2inv.d <- (1/la[2])*G2inv.n
      ed2 <- sum(dZtNZ*(G2inv.d*G^2))
      ed2 <- ifelse(ed2 == 0, 1e-50,ed2)
      tau2 <- sum(b.random^2*G2inv.n)/ed2
      tau2 <- ifelse(tau2 == 0, 1e-50, tau2)

      # Tau 3
      G3inv.d <- (1/la[3])*G3inv.n
      ed3 <- sum(dZtNZ*(G3inv.d*G^2))
      ed3 <- ifelse(ed3 == 0, 1e-50,ed3)
      tau3 <- sum(b.random^2*G3inv.n)/ed3
      tau3 <- ifelse(tau3 == 0, 1e-50, tau3)

      # New variance components and convergence check
      lanew <- c(tau1, tau2, tau3)
      #dla <- mean(abs(la - lanew))
      dla <- mean(abs((la - lanew)/lanew))
      la <- lanew

      if (trace) {
        cat(sprintf("%1$3d %2$10.6f", it, dla))
        cat(sprintf("%8.3f", c(ed1, ed2, ed3)), '\n')
      }
      if (dla < (thr[2])) break
    }

    if (trace) {
      end <- proc.time()[3]
      cat("Timings:\nSAP", (end - start), "seconds\n")
    }

    eta.old <- eta
    eta <- sclm::Xtheta(X, b.fixed) + sclm::Ztheta(Z, b.random, np[-1])

    gamma <- c(efine*exp(eta))
    mu <- sclm::Xtheta(C, gamma)

    # Convergence criterion (for trend)
    tol <- sum((eta - eta.old)^2)/sum(eta^2)
    print(tol)
    if (tol < (thr[1])) break
  }

  # Deviance, Effective dimension, AIC, and BIC
  dev <- 2*sum(y*log(ifelse(y == 0, 1, y/mu)) - (y - mu))
  ed <- prod(pord) + sum(c(ed1, ed2, ed3))
  aic <- dev + 2 * ed
  bic <- dev + log(length(y)) * ed

  opt.Ginv <- c((1/la[3])*g3u, (1/la[2])*g2u, (1/la[1])*g1u, (1/la[2])*g22b + (1/la[3])*g32b, (1/la[1])*g12b + (1/la[3])*g31b, (1/la[1])*g11b + (1/la[2])*g21b, (1/la[1])*g1t + (1/la[2])*g2t + (1/la[3])*g3t)

  end.all <- proc.time()[3]
  comp.time <- end.all - start.all

  cat("Computing time for all process:", comp.time, "seconds\n")

  output <- list(ndx = ndx, bdeg = bdeg, pord = pord, MM1 = MM1, MM2 = MM2, MM3 = MM3, var.comp = la, edf = c(ed1, ed2, ed3), b.fixed = b.fixed, b.random = b.random, eta = eta, gamma = gamma, mu = mu, maxit = i, computing.time = comp.time, dev = dev, ed = ed, aic = aic, bic = bic, matlist = list(X = X, Z = Z, C = C), Ginv = opt.Ginv, y = y, efine = efine)
  return(output)
}
