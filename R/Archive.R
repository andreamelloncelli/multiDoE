Archive <- function(dim, capacity) {
  ar <- list()
  if (nargs() < 2) {
    capacity <- 100
  }
  ar$nsols <- 0
  ar$dim <- dim
  ar$scores <- matrix(0, capacity, dim)
  ar$solutions <- vector("list", capacity)
  return(ar)
}

Resize <- function(ar) {
  n <- dim(ar$scores)[1]
  m <- dim(ar$scores)[2]
  ar$scores <- rbind(ar$scores, matrix(0, n, m))
  ar$solutions <- c(ar$solutions, vector("list", n))
  return(ar)
}

Add <- function(ar, sol, score) {
  # check for capacity
  if (length(ar$solutions) == ar$nsols) {
    ar <- Resize(ar)
  }
  ar$nsols <- ar$nsols + 1

  # add solution and scores
  ar$solutions[[ar$nsols]] <- sol
  ar$scores[ar$nsols, ] <- score
  return(ar)
}

# problema con repmat /kronecker
# prende in input una matrice
FixRepmat <- function(data, num_rows, num_cols){
  if (length(data) == 0) {
    return(matrix(0, 0, dim(data)[2] * num_cols))
  } else {
    return(repmat(data, num_rows, num_cols))
  }
}

RemoveDominated <- function(ar) {
  # select dominated solutions (and empty lines)
  toRemove <- matrix(0, ar$nsols, 1)

  for (i in 1:ar$nsols) {
    toRemove <- toRemove |
    apply((ar$scores >= FixRepmat(ar$scores[i, ], ar$nsols, 1)) &
      (ar$scores != FixRepmat(ar$scores[i, ], ar$nsols, 1)), 1, function(x) all(x))
  }
  # remove selected
  ar$solutions <- ar$solutions[ ! toRemove]
  ar$scores <- ar$scores[ ! toRemove, ]
  ar$nsols <- length(ar$solutions)
  return(ar)
}

RemoveDuplicates <- function(ar) {
  i <- ! duplicated(ar$scores)
  ar$scores <- ar$scores[i, ]
  ar$solutions <- ar$solutions[i]
  ar$nsols <- length(ar$solutions)
  return(ar)
}

Trim <- function(ar) {
  toRemove <- ! apply(ar$scores, 1, function(x) all(x == 0))
  ar$scores <- ar$scores[toRemove, ]
  ar$solutions <- ar$solutions[toRemove]
  ar$nsols <- length(ar$solutions)
  return(ar)
}





