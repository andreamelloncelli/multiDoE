#' @export
# function MSOpt ---------------------------

# 'units' e 'levels' nel codice matlab sono vettori (qui liste)

MSOpt = function(facts, units, levels, etas, criteria, model) {

  mso <- list()
  mso$facts <- facts
  mso$nfacts <- length(unlist(facts))
  mso$nstrat <- length(facts)
  mso$units <- units
  mso$runs <- prod(unlist(units))
  mso$etas <- etas
  mso$avlev <- as.list(rep(NA, mso$nfact))

  if (length(levels) == 1) {                                # is.scalar()
    mso$levs <- rep(1, mso$nfacts) * levels
    for (i in 1:mso$nfacts) {
      mso$avlev[[i]] <- (2 * 0:(levels - 1) / (levels - 1)) - 1
    }
  } else {
    mso$levs <- levels
    for (i in 1:mso$nfacts) {
      mso$avlev[[i]] <- (2 * 0:(levels[[i]] - 1) / (levels[[i]] - 1)) - 1
    }
  }

  V <- diag(mso$runs)
  for (i in 1:(mso$nstrat - 1)) {
    ones_shape <- prod(unlist(units)[(i + 1):length(units)])
    V <- V + etas[[1]] * kronecker(diag(prod(unlist(units)[1:i])),
                                   matrix(1, ones_shape, ones_shape))
  }

  mso$Vinv <- t(solve(V))
  mso$model <- model
  mso$crit <- criteria
  mso$ncrit <- length(criteria)

  if ("I" %in% criteria) {
    k <- mso$nfacts
    k2 <- k * (k - 1) / 2
    switch (model,
            "main" = {
              mso$M <- rbind(
                cbind(1, t(integer(k))),
                cbind(integer(k), diag(k) / 3)
              )
            },
            "interaction" = {
              mso$M <- rbind(
                cbind(1, t(integer(k)), t(integer(k2))),
                cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
                cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
              )
            },
            "quadratic" = {
              mso$M <- rbind(
                cbind(1, t(integer(k)), t(rep(1, k)) / 3, t(integer(k2))),
                cbind(integer(k), diag(k) / 3, matrix(0, k, k), matrix(0, k, k2)),
                cbind(rep(1, k) / 3, matrix(0, k, k),
                      (4 * diag(k) + 5 * matrix(1, k, k)) / 45, matrix(0, k, k2)),
                cbind(integer(k2), matrix(0, k2, k), matrix(0, k2, k), diag(k2) / 9)
              )
            },
            stop("Model type not valid")
    )
  }

  if ("Id" %in% criteria) {
    k <- mso$nfacts
    k2 <- k * (k - 1) / 2
    switch(model,
           "main" = {
             mso$M0 <- rbind(
               cbind(1, t(integer(k))),
               cbind(integer(k), diag(k) / 3)
             )
           },
           "interaction" = {
             mso$M0 <- rbind(
               cbind(1, t(integer(k)), t(integer(k2))),
               cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
               cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
             )
           },
           "quadratic" = {
             mso$M0 <- rbind(
               cbind(1, t(integer(k)), t(rep(1, k)) / 3, t(integer(k2))),
               cbind(integer(k), diag(k) / 3, matrix(0, k, k), matrix(0, k, k2)),
               cbind(rep(1, k) / 3, matrix(0, k, k),
                     (4 * diag(k) + 5 * matrix(1, k, k)) / 45, matrix(0, k, k2)),
               cbind(integer(k2), matrix(0, k2, k), matrix(0, k2, k), diag(k2) / 9)
             )
           },
           stop("Model type not valid")
    )
    mso$M0[, 1] <- 0
    mso$M0[1, ] <- 0
  }

  if ("As" %in% criteria) {
    switch(model,
           "main" = {
             w <- t(rep(1, mso$nfacts))
           },
           "interaction" = {
             w <- t(rep(1, mso$nfacts + mso$nfacts * (mso$nfacts - 1) / 2))
           },
           "quadratic" = {
             w <- c(
               rep(1, mso$nfacts),
               rep(1, mso$nfacts) / 4,
               rep(1, mso$nfacts * (mso$nfacts - 1) / 2)
             )
           },
           stop("Model type not valid")
    )
    a <- length(w / sum(w))
    mso$W <- c(w / sum(w)) * diag(a)
  }
  return(mso)
}

#' @export
# function colprod ---------------------------

colprod <- function(X) {

  n_row <- nrow(X)
  n_col <- ncol(X)
  out <- matrix(NA, n_row, n_col * (n_col - 1) / 2)

  k <- 1
  for (i in 1:(n_col - 1)) {
    for (j in (i + 1):n_col) {
      out[, k] <- X[, i] * X[, j]
      k <- k + 1
    }
  }
  return(out)
}

#' @export
# function Score ---------------------------

Score <- function(mso, settings) {

  switch(mso$model,
         "main" = {
           X <- cbind(rep(1, mso$runs), settings)
         },
         "interaction" = {
           X <- cbind(rep(1, mso$runs), settings, colprod(settings))
         },
         "quadratic" = {
           X <- cbind(rep(1, mso$runs), settings, settings ^ 2, colprod(settings))
         },
  )

  B <- t(X) %*% mso$Vinv %*% X
  determ <- det(B)
  scores <- as.vector(matrix(Inf, length(mso$crit)))

  if (rcond(B) > 1e-5 & determ > 0) {

    ind <- mso$crit == "D"                              # true/false
    if (any(ind)) {
      scores[ind] <- 1 / determ ^ (1 / dim(X)[2])
    }

    if (any(c("I", "Id", "Ds", "A", "As") %in% mso$crit)) {
      Binv <- solve(B)
    }

    ind <- mso$crit == "I"
    if (any(ind)) {
      scores[ind] <- sum(diag(Binv %*% mso$M))
    }

    ind <- mso$crit == "Id"
    if (any(ind)) {
      scores[ind] <- sum(diag(Binv %*% mso$M0))
    }

    ind <- mso$crit == "A"
    if (any(ind)) {
      scores[ind] <- sum(diag(Binv)) / dim(X)[2]
    }

    rws <- dim(Binv)[1]
    cls <- dim(Binv)[2]

    ind <- mso$crit == "Ds"
    if (any(ind)) {
      scores[ind] <- (det(Binv[2:rws, 2:cls])) ^ ( 1 / (dim(X)[2] - 1))
    }

    ind <- mso$crit == "As"
    if (any(ind)) {
      scores[ind] <- sum(diag(mso$W %*% Binv[2:rws, 2:cls]))
    }

  }
  return(scores)
}


