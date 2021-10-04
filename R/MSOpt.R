#' MSOpt
#'
#' @description The \code{MSOpt} function creates a list object containing
#' the main information on the experiment settings and the optimization
#' criteria to be considered for the optimal design construction.
#' According to the declared criteria, it also provides the basic matrices for
#' their implementation. \code{MSOpt} returns input objects of the
#' \code{\link[multiDoE]{Score}} and \code{\link[multiDoE]{MSSearch}} functions
#' of the multiDoE package.
#'
#' @usage MSOpt(facts, units, levels, etas, criteria, model)
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
#' @param criteria A list specifying the criteria to be optimized. It can
#' contain any combination of:
#' \itemize{
#'   \item{"I" : I-optimality}
#'   \item{"Id" : Id-optimality}
#'   \item{"D" : D-optimality}
#'   \item{"A" : Ds-optimality}
#'   \item{"Ds" : A-optimality}
#'   \item{"As" : As-optimality}
#' }
#' More detailed information on the available criteria is given under 'Details'.
#'
#' @param model A string which indicates the type of model, among "main",
#' "interaction" and "quadratic".
#'
#' @details A little notation is introduced to present the criteria
#' that can be used in the multi-objective approach of the multiDoE package. \cr
#' For an experiment with \eqn{N} runs and \eqn{s} strata, with stratum \eqn{i}
#' having \eqn{n_i} units within each unit at previous stratum \eqn{(i-1)} and
#' stratum 0 being defined as the entire experiment \eqn{(n_0 = 1)}, the
#' general form of the model can be written as:
#' \deqn{y = X\beta + \sum\limits_{i = 1}^{s} Z_i\varepsilon_i}
#' where \eqn{y} is a \eqn{N}-dimensional vector of responses
#' (\eqn{N = \prod_{j = 1}^{s}n_j}), \eqn{X} is an \eqn{N \times p} model matrix,
#' \eqn{\beta} is a \eqn{p}-dimensional vector containing the \eqn{p} fixed model
#' parameters, \eqn{Z_i} is an \eqn{N \times b_i} indicator matrix of zero and
#' ones for the units in stratum \eqn{i} and \eqn{b_i = \prod_{j = 1}^{i}n_j}.
#' Finally, the vector \eqn{\varepsilon_i \sim N(0,\sigma_i^2I_{b_i})} is a
#' \eqn{b_i}-dimensional vector containing the random effects, which are all
#' incorrelated. The variance components \eqn{\sigma^{2}_{i} (i = 1, \dots, s)}
#' have to be estimated and this is usually done by using the REML method.
#' The best linear unbiased estimator for the parameter vector \eqn{\beta} is
#' the generalized least square estimator:
#' \deqn{\hat{\beta}_{\emph{GLS}} = (X'*V^{-1}*X)^{-1}*X'*V^{-1}*y}.
#' This estimator has variance-covariance matrix:
#' \deqn{Var(\hat{\beta}_{\emph{GLS}) = \sigma^{2}(X'V^{-1}X)^{-1}}},
#' where \eqn{V = \sum\limits_{i = 1}^{s}\eta_i Z_i'Zi},
#' \eqn{\eta_i = \frac{\sigma_i^{2}}{\sigma^{2}}} and \eqn{\sigma^{2} =
#' \sigma^{2}_{s}}.
#' The variance components have to be estimated
#'
#' \itemize{
#'   \item \strong{\emph{I}-optimality.}
#'   \item \strong{\emph{Id}-optimality.}
#'   \item \strong{\emph{A}-optimality.}
#'   \item \strong{\emph{D}-optimality.}
#'   \item \strong{\emph{As}-optimality.}
#'   \item \strong{\emph{Ds}-optimality.}
#' }
#'
#' @return \code{MSOpt} returns a list containing the following components:
#' \item{facts}{The argument \code{facts}.}
#' \item{nfacts}{An integer. The number of expermental factors.}
#' \item{nstrat}{An integer. The number of strata.}
#' \item{units}{The argument \code{units}.}
#' \item{runs}{An integer. The number of runs.}
#' \item{etas}{The argument \code{etas}.}
#' \item{avlev}{A list showing the available levels for each experimental factor.}
#' \item{levs}{A vector showing the number of available levels for each experimental factor.}
#' \item{Vinv}{The inverse of the variance-covariance matrix of the responses.}
#' \item{model}{The argument \code{model}.}
#' \item{crit}{The argument \code{criteria}.}
#' \item{ncrit}{An integer. The number of criteria.}
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
#' The \code{Score} function returns the optimization criteria values for the
#' given \code{MSOpt} object and design matrix.
#'
#' @param msopt a list, output of \code{MSOpt} function.
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


