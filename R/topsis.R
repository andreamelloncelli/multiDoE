
Ldist <- function(x, y, w, p) {
  ldist <- ((w ^ p) %*% (abs(x - y) ^ p)) ^ (1 / p)
  return(ldist)
}

topsisOpt <- function(paretoFront, w, p) {

  if (sum(w) != 1) {
    stop("Weight vector invalid.")

  } else {
    matrice <- paretoFront$scores

    # 1. normalized score matrix
    rw <- dim(matrice)[1]
    cl <- dim(matrice)[2]
    normPF <- matrix(NA, rw, cl)
    normPF <- apply(matrice, 2, function(x) sapply(x, function(y) y/sum(x)))

    # 2. ideal and nadir point
    ideal <- apply(normPF, 2, min)
    nadir <- apply(normPF, 2, max)

    # 3. calculation of distances
    ideal_dist <- apply(normPF, 1, function(x) Ldist(x, ideal, w, p))
    nadir_dist <- apply(normPF, 1, function(x) Ldist(x, nadir, w, p))

    # 4. ranking
    S <- vector(length = rw)
    for (i in 1:rw) {
      S[i] <- nadir_dist[i] / (ideal_dist[i] + nadir_dist[i])
    }
    ranking <- data.frame("paretoFront_pos" = c(1:rw), "S(x)" = S)
    ranking2 <- ranking[order(ranking$S.x., decreasing = T),]
    bestScore <- paretoFront$scores[ranking2[1, 1],]
    bestSol <- paretoFront$solutions[[ranking2[1, 1]]]
  }
  return(list("ranking" = ranking2, "bestScore" = bestScore, "bestSol" = bestSol))
}

#
# paretoFront <- tpls$megaAR
# prova <- topsisOpt(paretoFront, w = c(1/5, 1/5, 3/5), p = 2)
# esse <- ranking[,2]
# plot_ly(x=paretoFront$scores[,1],
#         y=paretoFront$scores[,2],
#         z=paretoFront$scores[,3],
#         type="scatter3d", mode="markers", color = esse)

