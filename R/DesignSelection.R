#' optMultiCrit
#'
#' @description The \code{optMultiCrit} function suggests an objective criterion
#' for the selection of the best experimental design among all Pareto front solutions.
#' The selection is based on minimizing the euclidean distance in the criteria space
#' between all the Pareto front designs and an approximate utopian point. By default
#' the utopian point coordinates are the minimum value reached by every criteria
#' during an optimization procedure (\code{\link[multiDoE]{runTPLS}}); otherwise
#' it can be set to a specific value by the user.
#'
#' @param ar A list as the \code{megaAR} returned by the \code{runTPLS} function.
#' @param ... optional argument (see below).
#'
#' @details Additional arguments can be specified as follows:
#' \itemize{
#' \item \code{myUtopianPoint}: A vector containing the utopian point coordinates.
#' }
#'
#' @return The \code{optMultiCrit} function returns a list whose elements are:
#' \itemize{
#' \item \code{solution}: The selected optimal design matrix.
#' \item \code{score}: A vector containing the criteria scores for \code{solution}.
#' }
#
#' @export



optMultiCrit <- function(ar, ...) {

  varargin <- list(...)

  if (nargs() == 1) {
    bestPoint <- apply(ar$scores, 2, min)
  } else {
    bestPoint <- varargin[[1]]
  }

  d <- c()
  for (i in 1:ar$nsols) {
    d[i] <- dist(rbind(ar$scores[i, ], bestPoint))
  }

  ind <- which(d == min(d))
  return(list("solution" = ar$solutions[ind], "scores" = ar$score[ind, ]))
}

#' optSingleCrit
#'
#' @description The \code{optSingleCrit} function selects the Pareto front designs
#' that optimizes the individually considered criteria.
#'
#' @param ar A list as the \code{megaAR} returned by the \code{runTPLS} function.
#' @param criteria The criteria vector as the input of \code{runTPLS} function.
#'
#' @return A list whose \eqn{i}-th element corresponds to the solution that optimizes
#' the \eqn{i}-th criterion in \code{criteria}. The solution is a list of two elements:
#' \itemize{
#' \item \code{score}: A vector containing the scores for every element in \code{criteria}.
#' \item \code{solution}: The design matrix.
#' }
#' @export
#'
optSingleCrit <- function(ar) {

  nCrit <- dim(ar$scores)[2]
  best <- vector("list", nCrit)
  index <- apply(ar$scores, 2, which.min)

  for (i in 1:nCrit) {
    best[[i]] <- list(ar$scores[index[i],i], ar$solutions[[index[i]]])
  }

  names(best) <- colnames(ar$scores)
  return(best)
}


#' plotPareto
#'
#' @description The \code{plotPareto} function returns a graphical representation
#' (at most 3D) of the Pareto front.
#'
#' @usage plotPareto(ar, x, y)
#'
#' plotPareto(ar, x, y, z, mode = True)
#'
#' @param ar A list as the \code{megaAR} returned by the \code{runTPLS} function.
#' @param x The criterion on the x axis. It can be one of the following: \code{"I",
#' "Id", "D", "Ds", "A"} and \code{"As"}.
#' @param y The criterion on the y axis. It can be one of the following: \code{"I",
#' "Id", "D", "Ds", "A"} and \code{"As"}.
#' @param ... optional argument (see below).
#'
#' @details If three criteria are considered:
#' \itemize{
#' \item \code{z}: The criterion on the z axis. It can be one of the following:
#' \code{"I", "Id", "D", "Ds", "A"} and \code{"As"}.
#' \item \code{mode}: When \code{mode=True} the function returns a 3D interactive
#' chart. When \code{mode=False} it returns a 2D chart in which the \code{z} criteria
#' values are represented by a color scale.
#' }
#'
#' @return The Pareto front chart.
#' @export

plotPareto <- function(ar) {

  if (ar$nsols == 1) {
    nCrit <- length(ar$scores)
    df <- as.data.frame(t(ar$scores))
  } else {
    nCrit <- dim(ar$scores)[2]
    df <- as.data.frame(ar$scores)
  }

  if (nCrit == 2) {
    ggplot(df, aes_string(x = colnames(df)[1],
                    y = colnames(df)[2])) + geom_point()
  } else if (nCrit == 3) {
    print(df)
    scatterplot3d(df[, 1], df[, 2], df[,3])
    # rgl::scatter3d(df[, 1], df[, 2], df[,3])

  } else {
    stop("Number of criteria not valid")
  }
}


