#' MSOpt
#'
#' \code{MSOpt} initializes the framework for the implementation of several
#' optimization criteria for a multistratum design.
#'
#' @param facts a list representing the distribution of factors over strata.
#' @param units a list containing the number of units in each stratum.
#' @param levels a vector containing number of levels for each factor.
#' @param etas a list specifying ratios of error variance between subsequent
#' strata. The length must be equal to number of strata minus one.
#' @param criteria  the list of criteria to be optimized.
#' @param model a string which indicates the model type, among "main",
#' "interaction" and "quadratic".
#'
#'@details \code{criteria} can contain any combination of:
#' \itemize{\item "I" - I-optimality
#' \item "Id" - Id-optimality
#' \item "D" - D-optimality
#' \item "Ds" - Ds-optimality
#' \item "A" - A-optimality
#' \item "As" - As-optimality}
#'
#' See Sambo et al. (2016) for information on the available criteria.
#'
#' @return \code{MSOpt} returns a list, whose elements are:
#' \itemize{\item \code{facts} - a list representing the distribution of factors over strata.
#' \item \code{nfacts} - number of factors.
#' \item \code{nstrat} - number of strata.
#' \item \code{units} - a list containing the number of units in each stratum.
#' \item \code{runs} - number of runs.
#' \item \code{etas} - a list specifying ratios of error variance between subsequent strata.
#' \item \code{avlev} - a list showing the available levels for each factor.
#' \item \code{levs} - a vector showing the number of levels for each factor.
#' \item \code{Vinv} - the inverse of the covariance matrix of the responses.
#' \item \code{model} - a string indicating the model type.
#' \item \code{crit} - a vector containing the criteria to be optimize.
#' \item \code{ncrit} - number of criteria.
#' \item \code{M} - the matrix of moments of the cube. Only with I-optimality criteria.
#' \item \code{M0} - the matrix of moments of the cube. Only with Id-optimality criteria.
#' \item \code{W} - the diagonal matrix of weights. Only with As-optimality criteria.
#' }
#'
#' @export

MSOpt <- function(facts, units, levels, etas, criteria, model) {

  msopt <- list()

  msopt$facts <- facts
  msopt$nfacts <- length(unlist(facts))
  msopt$nstrat <- length(facts)
  msopt$units <- units
  msopt$runs <- prod(unlist(units))
  msopt$etas <- etas
  msopt$avlev <- as.list(rep(NA, msopt$nfact))

  if (length(levels) == 1) {
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
    if (i + 1 > length(units)) {
      ones_shape <- 1
    } else {
      ones_shape <- prod(unlist(units)[(i + 1):length(units)])
    }
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


# function colprod ---------------------------

colprod <- function(X) {

  if (is.matrix(X)) {
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
  } else if (is.vector(X) & length(X) > 1) {

    n <- length(X)
    out <- c()

    k <- 1
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        out[k] <- X[i] * X[j]
        k <- k + 1
      }
    }
  } else {
    out = NULL
  }
  return(out)
}


#' Score
#'
#' Scoring function for the given MSOpt list and design matrix. It contains the
#' implementation of all the available criteria.
#'
#' @param msopt a list, the output of MSOpt function.
#' @param settings the design matrix for which criteria scores are calculated.
#'
#' @return Vector of optimization criteria values.
#' @export

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
       scores[ind] <- round(1 / determ ^ (1 / dim(X)[2]), 10)
    }

    if (any(c("I", "Id", "Ds", "A", "As") %in% msopt$crit)) {
      Binv <- solve(B)
    }

    ind <- msopt$crit == "I"
    if (any(ind)) {
      # TODO remove round
      scores[ind] <- round(sum(diag(Binv %*% msopt$M)), 10)
    }

    ind <- msopt$crit == "Id"
    if (any(ind)) {
      # TODO remove round
      scores[ind] <- round(sum(diag(Binv %*% msopt$M0)), 10)
    }

    ind <- msopt$crit == "A"
    if (any(ind)) {
      # TODO remove round
      scores[ind] <- round(sum(diag(Binv)) / dim(X)[2], 10)
    }

    ind <- msopt$crit == "Ds"
    if (any(ind)) {
      rws <- dim(Binv)[1]
      cls <- dim(Binv)[2]
      # TODO remove round
      scores[ind] <- round((det(Binv[2:rws, 2:cls])) ^ ( 1 / (dim(X)[2] - 1)), 10)
    }

    ind <- msopt$crit == "As"
    if (any(ind)) {
      rws <- dim(Binv)[1]
      cls <- dim(Binv)[2]
      # TODO remove round
      scores[ind] <- round(sum(diag(msopt$W %*% Binv[2:rws, 2:cls])), 10)
    }
  }
  return(scores)
}


