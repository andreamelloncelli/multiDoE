
runTPLS <- function(facts, units, criteria, model, iters, ...) {

  optional = list(...)

  ar <- vector(mode = "list", iters)
  stats <- vector(mode = "list", iters)

  for (i in 1:iters) {
    tpls <- TPLSearch(facts, units, criteria, model, optional)
    ar[[i]] <- tpls$ar
    stats[[i]] <- tpls$stats
  }

  lCrit <- length(criteria)
  megaAR <- Archive(lCrit, iters * (restarts - lCrit  * (restInit - 1)))

  for (i in 1:iters) {
    for (j in 1: ar[[i]][["nsols"]]) {
      megaAR <- Add(ar[[i]][["solutions"]][[j]], ar[[i]][["scores"]][j, ])
    }
  }

  megaAR <- RemoveDuplicates(megaAR)
  megaAR <- RemoveDominated(megaAR)

  return(megaAR)
}
