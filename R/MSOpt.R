#' MSOpt
#'
#' @description The \code{MSOpt} function creates a list object containing
#' the main information on the experiment settings and the optimization
#' criteria to be considered. According to the declared criteria, it also
#' provides the basic matrices for their implementation. \code{MSOpt} returns
#' input objects of the \code{\link[multiDoE]{Score}} and
#' \code{\link[multiDoE]{MSSearch}} functions of the multiDoE package.
#'
#' @param facts A list of vectors representing the distribution of factors
#' across strata. Each item in the list represents a stratum and the first item
#' is the highest stratum of the multi-stratum structure of the experiment.
#' Within the vectors, experimental factors are indicated by progressive integer
#' from 1 to the total number of experimental factors, starting from the highest
#' strata. Blocking factors are denoted by empty vectors.
#'
#' @param units A list containing the number of units in each stratum. For
#' stratum i, the number of experimental units within each unit of the previous
#' stratum (i-1) is indicated. \code{length(units)} must be equal to
#' \code{length(facts)}.
#'
#' @param levels A vector containing the number of available levels for each
#' experimental factor (blocking factors are excluded). If all the experimental
#' factors share the number of levels, one integer is sufficient.
#'
#' @param etas A list specifying the ratios of error variance between subsequent
#' strata. It follows that \code{length(etas)} must be equal to
#' \code{length(facts) - 1}.
#'
#' @param criteria A list containing the criteria to be optimized. It can
#' contain any combination of:
#' \itemize{
#'   \item{"I" : I-optimality}
#'   \item{"Id" : Id-optimality}
#'   \item{"D" : D-optimality}
#'   \item{"Ds" : Ds-optimality}
#'   \item{"A" : A-optimality}
#'   \item{"As" : As-optimality}
#' }
#' See the \strong{Details} section for more detailed information on the available
#' criteria. \cr
#' The details of model specification are given under 'Details'.
#'
#' @param model A string which indicates the type of model, among "main",
#' "interaction" and "quadratic".
#'
#' @details
#' In order to... \cr
#' The general form of the model:
#' \deqn{}
#'
#' Generalized least square estimator:
#' \deqn{\hat{\beta}_{\emph{GLS}} = (X'*V^{-1}*X)^{-1}*X'*V^{-1}*y}
#'
#' \itemize{
#'   \item \strong{\emph{D}-optimality.} The
#'
#'   \item {"Ds"} {Ds-optimality}
#'   \item {"A"} {A-optimality}
#'   \item {"As"} {As-optimality}
#' }
#'
#' @return \code{MSOpt} returns a list containing the following components:
#' \item{facts}{The argument \code{facts}.}
#' \item{nfacts}{An integer, the number of expermental factors.}
#' \item{nstrat}{An integer, the number of strata.}
#' \item{units}{The argument \code{units}.}
#' \item{runs}{An integer, the number of runs.}
#' \item{etas}{The argument \code{etas}.}
#' \item{avlev}{A list showing the available levels for each experimental factor.}
#' \item{levs}{A vector showing the number of levels for each experimental factor.}
#' \item{Vinv}{The inverse of the variance-covariance matrix of the responses.}
#' \item{model}{The argument \code{model}.}
#' \item{crit}{The argument \code{criteria}.}
#' \item{ncrit}{An integer, the number of criteria.}
#' \item{M}{The matrix of moments of the cube. Only with \emph{I-optimality} criteria.}
#' \item{M0}{The matrix of moments of the cube. Only with \emph{Id-optimality} criteria.}
#' \item{W}{The diagonal matrix of weights. Only with \emph{As-optimality} criteria.}
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
    msopt$levs <- rep(levels, msopt$nfacts)
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
       scores[ind] <- 1 / determ ^ (1 / dim(X)[2])
    }

    if (any(c("I", "Id", "Ds", "A", "As") %in% msopt$crit)) {
      Binv <- solve(B)
    }

    ind <- msopt$crit == "I"
    if (any(ind)) {
      scores[ind] <- sum(diag(Binv %*% msopt$M))
    }

    ind <- msopt$crit == "Id"
    if (any(ind)) {
      scores[ind] <- sum(diag(Binv %*% msopt$M0))
    }

    ind <- msopt$crit == "A"
    if (any(ind)) {
      scores[ind] <- sum(diag(Binv)) / dim(X)[2]
    }

    ind <- msopt$crit == "Ds"
    if (any(ind)) {
      rws <- dim(Binv)[1]
      cls <- dim(Binv)[2]
      scores[ind] <- (det(Binv[2:rws, 2:cls])) ^ ( 1 / (dim(X)[2] - 1))
    }

    ind <- msopt$crit == "As"
    if (any(ind)) {
      rws <- dim(Binv)[1]
      cls <- dim(Binv)[2]
      scores[ind] <- sum(diag(msopt$W %*% Binv[2:rws, 2:cls]))
    }
  }
  return(scores)
}


