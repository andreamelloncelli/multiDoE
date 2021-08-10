
# bestMultiCrit: Set disegni ottimi, ottimizzazione multicriterio
# Funzione che prende in input il risultato della funzione TPLSearch (un oggetto
# di classe Archive) che contiene le soluzioni appartenente al fronte di Pareto
# e le coordinate del punto utopico.
# La funzione restituisce una lista che ha come primo elemento una matrice che contiene
# nelle righe gli scores delle soluzioni ottime e come secondo elemento la lista delle
# soluzioni (disegni ottimi)

bestMultiCrit <- function(ar, point) {
  d <- c()
  for (i in 1:ar$nsols) {
    d[i] <- dist(rbind(ar$scores[i, ], point))
  }
  ind <- which(d == min(d))
  return(list("scores" = ar$score[ind, ], "solution" = ar$solutions[ind]))
}

# esempio: y è oggetto Archive
# y = list()
# y$nsols = 4
# y$scores = matrix(c(1,1,1,1,2,2,2,2,3,3,3,3, 1,1,1,1), nrow = 4, ncol = 4, byrow = T)
# y$solutions = list(1,2,3)
# bestDesignMC(y, c(1,1,1,1))


bestSingleCrit <- function(ar) {
  nCrit <- dim(ar$scores)[2]
  best <- vector("list", nCrit)
  names(best) <- colnames(ar$scores)
  index <- apply(ar$scores, 2, which.min)
  for (i in 1:nCrit) {
    best[[i]] <- list(ar$scores[index[i],i], ar$solutions[[index[i]]])
  }
  return(best)
}
#
# plotPareto <- function(ar) {
#   nCrit <- dim(ar$scores)[2]
#   if (nCrit == 2) {
#     df <- as.data.frame(ar$scores)
#     ggplot(data = df, aes(x = colnames(ar$scores)[1],
#                           y = colnames(ar$scores)[2]) + geom_point()
#   } else {
#
#   }
# }
