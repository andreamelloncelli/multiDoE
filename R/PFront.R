library(pracma)

PFront <- function(arch) {
  pf = list()

  pf$arch <- arch
  pf$ptrs <- c()    # NOTA: era nelle properties
  pf$gaps <- c()    # NOTA: era nelle properties
  pf$scmax <- matrix(0, 1, arch$dim)
  pf$scmin <- matrix(Inf, 1, arch$dim)

  # add the last arch.dim entries to the Pareto Front
  for (i in 0:(arch$dim - 1)) {
    Add_PF(pf, pf$arch$nsols - i)
  }
  UpdateMinMaxSc(pf)
  # UpdateGaps
  return(pf)
}

rowleq <- function(A, B) {
  argmin = which.min(A==B)
  return(A[argmin] <= B[argmin])
}

# Add è gia una funzione della classe Archive
# sostituisco con Add_PF
# arch$: nsols, dim, scores, solutions

Add_PF <- function(pf, solPtr) {

  if (IsWeakDominated(solPtr)) {
    return()
  }

  wDom <- GetWeakDominated(solPtr)
  pf$ptrs <- pf$ptrs(!wDom)

  # insert in the correct position
  if (length(pf$ptrs) == 0) {
    pf$ptrs <- solPtr
  } else if (rowleq(pf$arch$scores[pf$ptrs[length(pf$ptrs)], ], pf$arch$scores[solPtr, ])) {
    pf$ptrs <- rbind(pf$ptrs, solPtr)
  } else {
    for (i in 1:length(pf$ptrs)) {
      if (rowleq(pf$arch$scores[solPtr, ], pf$arch$scores[pf$ptrs[i], ])) {
        pf$ptrs = rbind(pf$ptrs[1:(i-1)], solPtr, pf$ptrs[i:length(pf$ptrs)])
        break
      }
    }
  }
  UpdateMinMaxSc(pf)
}

AddNoNorm <- function(pf, solPtr) {

  if (IsWeakDominated(solPtr)) {
    return()
  }

  wDom <- GetWeakDominated(solPtr)
  pf$ptrs <- pf$ptrs(!wDom)

  # insert in the correct position
  if (length(pf$ptrs) == 0) {
    pf$ptrs <- solPtr
  } else if (rowleq(pf$arch$scores[pf$ptrs[length(pf$ptrs)], ], pf$arch$scores[solPtr, ])) {
    pf$ptrs <- rbind(pf$ptrs, solPtr)
  } else {
    for (i in 1:length(pf$ptrs)) {
      if (rowleq(pf$arch$scores[solPtr, ], pf$arch$scores[pf$ptrs[i], ])) {
        pf$ptrs = rbind(pf$ptrs[1:(i-1)], solPtr, pf$ptrs[i:length(pf$ptrs)])
        break
      }
    }
  }
}

# require(mco)

HyperVolume <- function(pf) {
  hv <- dominated_hypervolume(pf$Getnorm(pf$ptrs),
                              matrix(1, pf$arch$dim, 1) * 1.1)
  return(hv)
}

IsWeakDominated <- function(pf, solPtr) {
  flag <- apply(apply(repmat(pf$arch$scores[solPtr, ], length(pf$ptrs), 1) >=
    pf$arch$scores[pf$ptrs, ], 1, all), 2, any)
  return(flag)
}

GetWeakDominated <- function(pf, solPtr) {
  wDom <- apply(pf$arch$scores[pf$ptrs, ] >=
    repmat(pf$arch$scores[solPtr, ], length(pf$ptrs), 1), 1, all)
  return(wDom)
}

UpdateMinMaxSc <- function(pf) {
  pf$scmax <- apply(pf$arch$scores[pf$ptrs, ], 2, max)
  pf$scmin <- apply(pf$arch$scores[pf$ptrs, ], 2, min)
}

SetMinMaxSc <- function(pf, scmax, scmin) {
  pf$scmax <- scmax
  pf$scmin <- scmin
}

GetNorm <- function(pf, solInd) {
  cnorm <- (pf$arch$scored[solInd, ] - repmat(pf$scmin, length(solInd), 1)) /
    repmat(pf$scmax - pf$scmin, length(solInd), 1)
  return(cnorm)
}
