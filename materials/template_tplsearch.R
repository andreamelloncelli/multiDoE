options(digits = 10)
# template TPLSearch
facts <- list(1:2, 3:4)
units <- list(12, 4)
levels <- 3
etas <- list(1, 1, 1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "quadratic"

lCrit <- length(criteria)
iters <- 10 * lCrit
restarts <- 100
restInit <- 2

ar <- vector(mode = "list", iters)
stats <- vector(mode = "list", iters)

#tStart <- tic

for (i in 1:iters) {
  print(i)
  tpls <- TPLSearch(facts, units, criteria, model, "Etas", etas,
                    "Levels", levels, "Restarts", restarts, "RestInit",
                    restInit, "RngSeed", i)
  ar[[i]] <- tpls$ar
  stats[[i]] <- tpls$stats
}






