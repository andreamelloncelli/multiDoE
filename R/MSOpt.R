#' MSOpt
#'
#' description
#'
#' @param facts a list of elements representing the distribution of factors over strata.
#'
#' @param units a list of numbers specifying the number of units in each stratum.
#' @param levels a list of numbers specifying the number of levels for each factor.
#' @param etas a list containing ratios of error variance between subsequent strata
#' @param criteria list of criteria to be optimized. It can contain any combination of:
#' \itemize{\item "I" - I-optimality
#'          \item "Id" - Id-optimality}
#'
#' \deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{%
#' p(x) = \lambda^x exp(-\lambda)/x!}
#' for \eqn{x = 0, 1, 2, \ldots}
#'
#' @param model type of model, among "main", "interaction" or "quadratic"
#'
#' @return MSOpt
#'
#' @export
#'
#'
MSOpt = function(facts, units, levels, etas, criteria, model) {

  msopt <- list()
  msopt$facts <- facts
  msopt$nfacts <- length(unlist(facts))
  msopt$nstrat <- length(facts)
  msopt$units <- units
  msopt$runs <- prod(unlist(units))
  msopt$etas <- etas
  msopt$avlev <- as.list(rep(NA, msopt$nfact))

  if (length(levels) == 1) {                      # NOTA: is.scalar()
    msopt$levs <- rep(1, msopt$nfacts) * levels
    for (i in 1:msopt$nfacts) {
      msopt$avlev[[i]] <- (2 * 0:(levels - 1) / (levels - 1)) - 1
    }
  } else {
    msopt$levs <- levels
    for (i in 1:msopt$nfacts) {
      msopt$avlev[[i]] <- (2 * 0:(levels[[i]] - 1) / (levels[[i]] - 1)) - 1
    }
  }

  V <- diag(msopt$runs)
  for (i in 1:(msopt$nstrat - 1)) {
    ones_shape <- prod(unlist(units)[(i + 1):length(units)])
    V <- V + etas[[1]] * kronecker(diag(prod(unlist(units)[1:i])),
                                   matrix(1, ones_shape, ones_shape))
  }

  msopt$Vinv <- t(solve(V))
  msopt$model <- model
  msopt$crit <- criteria
  msopt$ncrit <- length(criteria)

  if ("I" %in% criteria) {
    k <- msopt$nfacts
    k2 <- k * (k - 1) / 2
    switch (model,
           "main" = {
             msopt$M <- rbind(
               cbind(1, t(integer(k))),
               cbind(integer(k), diag(k) / 3)
               )
             },
           "interaction" = {
             msopt$M <- rbind(
               cbind(1, t(integer(k)), t(integer(k2))),
               cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
               cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
               )
             },
           "quadratic" = {
             msopt$M <- rbind(
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
    k <- msopt$nfacts
    k2 <- k * (k - 1) / 2
    switch(model,
           "main" = {
             msopt$M0 <- rbind(
               cbind(1, t(integer(k))),
               cbind(integer(k), diag(k) / 3)
               )
             },
           "interaction" = {
             msopt$M0 <- rbind(
               cbind(1, t(integer(k)), t(integer(k2))),
               cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
               cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
               )
             },
           "quadratic" = {
             msopt$M0 <- rbind(
               cbind(1, t(integer(k)), t(rep(1, k)) / 3, t(integer(k2))),
               cbind(integer(k), diag(k) / 3, matrix(0, k, k), matrix(0, k, k2)),
               cbind(rep(1, k) / 3, matrix(0, k, k),
                     (4 * diag(k) + 5 * matrix(1, k, k)) / 45, matrix(0, k, k2)),
               cbind(integer(k2), matrix(0, k2, k), matrix(0, k2, k), diag(k2) / 9)
               )
             },
           stop("Model type not valid")
    )
    msopt$M0[, 1] <- 0
    msopt$M0[1, ] <- 0
  }

  if ("As" %in% criteria) {
    switch(model,
           "main" = {
             w <- t(rep(1, msopt$nfacts))
             },
           "interaction" = {
             w <- t(rep(1, msopt$nfacts + msopt$nfacts * (msopt$nfacts - 1) / 2))
             },
           "quadratic" = {
             w <- c(
               rep(1, msopt$nfacts),
               rep(1, msopt$nfacts) / 4,
               rep(1, msopt$nfacts * (msopt$nfacts - 1) / 2)
               )
             },
           stop("Model type not valid")
    )
    a <- length(w / sum(w))
    msopt$W <- c(w / sum(w)) * diag(a)
  }
  return(msopt)
}



#' @export
#'
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


#' Score
#'
#' description
#'
#' @param msopt
#' @param settings
#'
#' @return
#' @export
#'
#' @examples
#'
Score <- function(msopt, settings) {

  switch(msopt$model,
         "main" = {
           X <- cbind(rep(1, msopt$runs), settings)
           },
         "interaction" = {
           X <- cbind(rep(1, msopt$runs), settings, colprod(settings))
           },
         "quadratic" = {
           X <- cbind(rep(1, msopt$runs), settings, settings ^ 2, colprod(settings))
           },
         )

  B <- t(X) %*% msopt$Vinv %*% X
  determ <- det(B)
  scores <- as.vector(matrix(Inf, length(msopt$crit)))


  if (rcond(B) > 1e-5 & determ > 0) {

    ind <- msopt$crit == "D"             # true/false
    if (any(ind)) {
       # TODO remove round
       scores[ind] <- round(1 / determ ^ (1 / dim(X)[2]),4)

    }

    if (any(c("I", "Id", "Ds", "A", "As") %in% msopt$crit)) {
      Binv <- solve(B)
    }

    ind <- msopt$crit == "I"
    if (any(ind)) {
      # TODO remove round
      scores[ind] <- round(sum(diag(Binv %*% msopt$M)), 4)
    }

    ind <- msopt$crit == "Id"
    if (any(ind)) {
      # TODO remove round
      scores[ind] <- round(sum(diag(Binv %*% msopt$M0)), 4)
    }

    ind <- msopt$crit == "A"
    if (any(ind)) {
      # TODO remove round
      scores[ind] <- round(sum(diag(Binv)) / dim(X)[2], 4)
    }

    ind <- msopt$crit == "Ds"
    if (any(ind)) {
      rws <- dim(Binv)[1]
      cls <- dim(Binv)[2]
      # TODO remove round
      scores[ind] <- round((det(Binv[2:rws, 2:cls])) ^ ( 1 / (dim(X)[2] - 1)), 4)
    }

    ind <- msopt$crit == "As"
    if (any(ind)) {
      # TODO remove round
      scores[ind] <- round(sum(diag(msopt$W %*% Binv[2:rws, 2:cls])), 4)
    }

  }
  print("funzione score")
  print(scores)
  return(scores)
}


