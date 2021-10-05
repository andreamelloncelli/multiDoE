#' runTPLS
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
#'
#' @param model A string which indicates the type of model, among "main",
#' "interaction" and "quadratic".
#'
#' @param iters An integer indicating the number of iterations of the MS-TPLS
#' algorithm.
#'
#' @param ... optional arguments (see the 'Details' section).
#'
#' @details Additional arguments can be specified as follows:
#' \item{'Rest arts', restarts}{A string and an integer, used in pair.}
#' \item{'Levels', levels}{Used in pair, the second parameter may be a
#' vector or an integer.}
#' \item{'Etas', etas}{Used in pair}
#'
#' \item{'RngSeed', rngSeed}{# ?}
#' \item{'StartMode', startMode}
#' \item{'AlphaMode', alphaMode}
#' \item{'RestInit', restInit}
#' \item{'Interactive', interact}
#'
#' @return
#' @export

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
