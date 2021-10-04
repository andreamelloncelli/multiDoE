#' MSSearch
#'
#' @description The \code{MSSearch} function implement the single-objective
#' local search component of the MS-TPLS algorithm.
#'
#'  is used to search for the optimal design, minimizing
#' for the following scalarization between the criteria:
#' \deqn{w = alpha * (Crit - CritTR) / CritSC}
#' where \emph{alpha} is the vector of the relative weights between the criteria;
#' \emph{Crit} is the vector of criteria values; \emph{CritTR} and \emph{CritSC}
#' are optional normalization factors.
#'
#' @usage MSSearch(msopt, alpha, "Start", sol, "Restarts", r, "Normalize",
#' c(CritTR, CritSC))
#'
#' @param msopt a list. The output of the MSOpt function.
#' @param alpha a vector of weights. The elements must add up to one.
#' @param ... optional arguments (see below)
#'
#' @details Additional arguments can be specified as follows:
#' \itemize{
#' \item \code{'Start', sol}: a string and a matrix, used in pair. They provide
#' a starting solution (\code{sol}) to the algorithm. By default initial solution
#' is randomly sampled.
#' \item \code{'Restarts', r }: a string and an integer, used in pair. They restart
#' the algorithm \code{r} times and finally the best solution is considered. By
#' default \eqn{r = 1}.
#' \item \code{'Normalize', c(CritTR, CritSC)}: a string and a vector. The second
#' is the vector the optional normalization factors. \code{CritTR} and \code{CritSC}
#' are vectors of length equal to the number of criteria, whose default elements
#' are 0 and 1 respectively.
#' }
#'
#' @return \code{MSSearch} returns a list, whose elements are:
#' \itemize{
#' \item \code{optsol} a design matrix. The best solution found.
#' \item \code{optscore} a vector containing the criteria scores.
#' \item \code{feval} an integer representing the number of score function
#' evaluations.
#' \item \code{trend} a vector of length \code{r} containing the minimum value of
#' \emph{w} for each iteration.
#' }
#'
#' @export


MSSearch <- function(msopt, alpha, ...) {
  varargin <- list(...)

  # default parameters
  restarts <- 1
  norms <- t(c(integer(msopt$ncrit), rep(1, msopt$ncrit)))
  random_start <- 1
  sol <- matrix(0, msopt$runs, msopt$nfacts)

  # optional parameters
  if (nargs() > 3) {
    for (i in seq(1, nargs() - 3, 2)) {
      switch (varargin[[i]],
              "Start" = {
                sol <- varargin[[i + 1]]
                random_start <- 0;
              },
              "Normalize" = {
                norms <- t(varargin[[i + 1]])
              },
              "Restarts" = {
                restarts = varargin[[i + 1]]
              }
      )
    }
  }

  CritTR <- norms[1:msopt$ncrit]
  CritSC <- norms[(msopt$ncrit + 1):length(norms)]

  # useful for working stratum by stratum
  totUnits <- t(numeric(msopt$ncrit))
  sizUnits <- totUnits

  for (s in 1:(msopt$nstrat-1)) {
    totUnits[s] <- prod(unlist(msopt$units[1:s]))
    sizUnits[s] <- prod(unlist(msopt$units[(s + 1):length(msopt$units)]))
  }
  totUnits[msopt$nstrat] <- prod(unlist(msopt$units[1:msopt$nstrat]))
  sizUnits[msopt$nstrat] <- 1

  optsc <- matrix(Inf, 1, msopt$ncrit)
  wopt <- Inf
  optsol <- matrix(0, msopt$runs, msopt$nfacts)
  feval <- 0
  trend <- numeric(restarts)

  for (t in 1:restarts) {
    if (random_start) {       # generate initial random solution

      # sample for each stratum
      for (s in 1:msopt$nstrat) {
        for (i in 1:totUnits[s]) {
          for (j in msopt$facts[[s]]) {
            sol[(sizUnits[s]*(i - 1) + 1):(sizUnits[s]*i), j] <-
              msopt$avlev[[j]][sample(1:msopt$levs[j], 1)]
          }
        }
      }
    } # if

    # score initial solution
    score <- Score(msopt, sol)
    wscore <- as.numeric(as.vector(alpha) %*% as.vector((score - CritTR) / CritSC))

    if (is.nan(wscore) | is.na(wscore)){
      wscore = Inf
    }
    feval <- feval + 1

    # Try to improve
    improvement <- 1
    sol2 <- sol

    while (improvement == 1) {
      improvement <- 0

      # Vary each factor in each stratum
      for (s in 1:msopt$nstrat) {
        for (i in 1: totUnits[s]) {
          for (j in msopt$facts[[s]]) {
            for (f in 1:length(msopt$avlev[[j]])) {

              if (sol[i * sizUnits[s], j] != msopt$avlev[[j]][f]) {
                sol2[(sizUnits[s]*(i - 1) + 1):(i * sizUnits[s]), j] <- msopt$avlev[[j]][f]
                score2 <- Score(msopt, sol2)

                wscore2 <- as.numeric(as.vector(alpha) %*% as.vector((score2 - CritTR) / CritSC))
                feval <- feval + 1
                #print(wscore2)
                if ((!is.nan(wscore2)) && (!is.na(wscore2)) && (wscore2 < wscore)){
                  sol <- sol2
                  wscore <- wscore2
                  improvement <- 1
                } else {
                  sol2 <- sol
                }
              }
            }
          }
        }
      }
    } # while

    if (wscore < wopt) {
      optsol <- sol
      optsc <- Score(msopt, sol)
      wopt <- wscore
    }

    trend[t] <- wopt
  } # for t
  return(list("optsol" = optsol, "optsc" = optsc, "feval" = feval, "trend" = trend))
}



