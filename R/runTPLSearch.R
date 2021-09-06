#' runTPLS
#'
#' @param facts a list representing the distribution of factors over strata.
#' @param units a list containing the number of units in each stratum.
#' @param criteria the list of criteria to be optimized.
#' @param model a string which indicates the model type, among "main",
#' "interaction" and "quadratic".
#' @param iters an integer indicating the number of iterations.
#' @param ... optional arguments (see below).
#'
#' @details Additional arguments can be specified as follows:
#' \itemize{
#' \item \code{'Restarts', restarts}
#' \item \code{'Levels', levels}
#' \item \code{'Etas', etas}
#' \item \code{'RngSeed', rngSeed}
#' \item \code{'StartMode', startMode}
#' \item \code{'AlphaMode', alphaMode}
#' \item \code{'RestInit', restInit}
#' \item \code{'Interactive', interact}
#' }
#'
#' @return
#' @export
#'

runTPLS <- function(facts, units, criteria, model, iters, ...) {

  varargin <- list(...)

  ar <- vector(mode = "list", iters)
  stats <- vector(mode = "list", iters)


  for (i in 1:iters) {
    print(i)
    varargin[which(varargin == "RngSeed") + 1] <- i
    tpls <- TPLSearch(facts, units, criteria, model, varargin)
    ar[[i]] <- tpls$ar
    stats[[i]] <- tpls$stats
  }

  lCrit <- length(criteria)
  megaAR <- Archive(lCrit, iters * (restarts - lCrit  * (restInit - 1)))

  for (i in 1:iters) {
    for (j in 1: ar[[i]][["nsols"]]) {
      megaAR <- Add(megaAR, ar[[i]][["solutions"]][[j]], ar[[i]][["scores"]][j, ])
    }
  }

  megaAR <- RemoveDuplicates(megaAR)
  megaAR <- RemoveDominated(megaAR)
# colnames(megaAR$scores) <- criteria

  return(list("ar" = ar, "stats" = stats, "megaAR" = megaAR))
}

# dare nomi agli elementi della lista in output
