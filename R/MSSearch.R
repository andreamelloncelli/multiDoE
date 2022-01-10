#' MSSearch
#'
#' @description The \code{MSSearch} function implements the extension of the
#' coordinate-exchange procedure proposed by Sambo, Borrotti, Mylona e Gilmour (2016)
#' called MS-Opt. This function can be used for the construction of optimal
#' multi-stratum experimental design considering one or more criteria (at most 6
#' criteria simultaneously).
#' The implemented algorithm seeks to minimize the following scalarization between criteria:
#' \deqn{f_W = \sum_{c \in C}{\alpha_cf_c(d; \eta)=\overline{\alpha} \cdot \overline{f}}, \quad \sum_{c \in C}\alpha_c = 1}{%
#' fW = \sum{c in C}{\alpha_c f_c(d;\eta)=\overline{\alpha} \cdot \overline{f}}, \sum{c in C}\alpha_c = 1,}
#' where \eqn{c} is the set of criteria to be minimized, \eqn{f_c} is the objective function
#' for the \eqn{c} criterion and \eqn{\overline{\alpha}} is the vector that controls
#' the relative weights of the objective functions. When there's a single criterion
#' of interest the function \code{MSSearch} corresponds to the single-objective local
#' search component of the MS-TPLS algorithm.
#'
#' @usage MSSearch(msopt, alpha, "Start", sol, "Restarts", r, "Normalize",
#' c(CritTR, CritSC )))
#'
#' @param msopt A list as returned by the \code{\link[multiDoE]{MSOpt}} function.
#' @param alpha A vector of weights, whose elements must add up to one.
#' \code{length(\eqn{\alpha})} must be equal to the number of criteria considered.
#' @param ... optional arguments (see below).
#'
#' @details Additional arguments can be specified as follows:
#' \itemize{
#' \item \code{'Start', sol}: a string and a matrix, used in pair. They provide
#' a starting solution (\code{sol}) or initial design to the algorithm. By
#' default the initial solution is randomly generated following the SampleDesign()
#' procedure described in Sambo, Borrotti, Mylona and Gilmour (2016).
#'
#' \item \code{'Restarts', r }: a string and an integer, used in pair. When \code{r=1},
#' the default value, the procedure implemented in \code{MSSearch} is a local search
#' algorithm optimizing the objective function \eqn{f_W} starting from one initial
#' design in the design space. This parameter allows to restart the algorithm \code{r} times.
#' If no initial design is passed, for each iteration, a different starting solution
#' is generated, letting the probability to find a global minimum be higher.
#' The design returned by the algorithm is the one that minimizes \eqn{f_W} across
#' all the iterations.
#'
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





