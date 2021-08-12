
runTPLS <- function(facts, units, criteria, model, iters, ...) {

  varargin <- list(...)

  ar <- vector(mode = "list", iters)
  stats <- vector(mode = "list", iters)

  for (i in 1:iters) {
    print(i)
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

  return(list(ar, stats, megaAR))
}

# dare nomi agli elementi della lista in output
