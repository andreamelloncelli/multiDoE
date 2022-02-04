
Ldist <- function(x, y, w, p) {
  ldist <- ((w ^ p) %*% (abs(x - y) ^ p)) ^ (1 / p)
  return(ldist)
}

#'
#'
#' @param matrice
#' @param w
#' @param p
#'
#' @return
#' @export
#'


#plot_ly(x=paretoFront$scores[,1],
#        y=paretoFront$scores[,2],
#        z=paretoFront$scores[,3],
#        type="scatter3d", mode="markers", color = esse)

