#' optMultiCrit
#'
#' \code{optMultiCrit} is used to select the best solution(s) from the optimal
#' experimental designs belonging to the Pareto front. The selection is based on
#' minimizing the Euclidean distance between the optimal solutions and an approximate
#' utopian point. The latter is calculated by default as the vector of the minimum
#' criteria values found during a multi-objective optimization algorithm (\code{runTPLS}
#' function in the \code{MultiDoE} package). Alternatively, a different point can
#' be chosen.
#'
#' @param ar a list. It is an archive object containing:
#' \itemize
#' \item an integer. The number of solutions.
#' \item an integer. The number of criteria considered.
#' \item a matrix. The score matrix associated with the solutions found.
#' \item a list. The list of solutions (design matrices).
#' @param ... optional argument (see below).
#'
#' @details
#'
#' @return \code{optMultiCrit} function returns a list. The first element is a list
#' of design matrices, the second the is the matrix of the respective criteria scores.
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
#' @param ar
#'
#' @return
#' @export
#'

optSingleCrit <- function(ar, criteria) {

  nCrit <- dim(ar$scores)[2]
  best <- vector("list", nCrit)
  names(best) <- criteria
  index <- apply(ar$scores, 2, which.min)

  for (i in 1:nCrit) {
    best[[i]] <- list(ar$scores[index[i],i], ar$solutions[[index[i]]])
  }
  return(best)
}


#' plotPareto
#'
#' @param ar
#'
#' @return
#' @export
#'

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


