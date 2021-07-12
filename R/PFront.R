# library(pracma)
library(mco)

PFront <- function(arch) {
  pf <- vector(mode = "list")

  pf$arch <- arch
  pf$ptrs <- c()
  pf$gaps <- list()
  pf$scmax <- matrix(0, 1, arch$dim)
  pf$scmin <- matrix(Inf, 1, arch$dim)

  # add the last arch.dim entries to the Pareto Front
  for (i in 0:(arch$dim - 1)) {
    pf <- Add_PF(pf, solPtr = pf$arch$nsols - i)
  }

  pf <- UpdateMinMaxSc(pf)
  # UpdateGaps
  return(pf)
}

# Add_PF perchè Add è gia una funzione della classe Archive

# arch$: nsols, dim, scores, solutions

# Add_PF e AddNoNorm (quasi uguali) possono essere unite in un'unica funzione con un
# argomento T/F opzionale ?

Add_PF <- function(pf, solPtr) {

  if (IsWeakDominated(pf, solPtr)) {
    return(pf)
  }

  # fino a qui ok

  wDom <- GetWeakDominated(pf, solPtr)
  print("wdommmmm")
  print(wDom)
  print("we")

  pf$ptrs <- pf$ptrs[!wDom]

  # se lo lasciassi fuori dal ciclo restituirebbe NULL
  print("pfptrs")
  print(pf$ptrs)

  # insert in the correct position
  if (length(pf$ptrs) == 0) {
    print("entra nel secondooooo if")
    print("solptr")
    print(solPtr)
    pf$ptrs <- solPtr
  } else if (rowleq(pf$arch$scores[pf$ptrs[length(pf$ptrs)], ],
                    pf$arch$scores[solPtr, ])) {
    print("entra nel terzo")
    pf$ptrs <- rbind(pf$ptrs, solPtr)
  } else {
    print("entra nel quarto")
    for (i in 1:length(pf$ptrs)) {
      if (rowleq(pf$arch$scores[solPtr, ], pf$arch$scores[pf$ptrs[i], ])) {
        pf$ptrs = rbind(pf$ptrs[1:(i-1)], solPtr, pf$ptrs[i:length(pf$ptrs)])
        break
      }
    }
  }

  pf <- UpdateMinMaxSc(pf)
  return(pf)
}

AddNoNorm <- function(pf, solPtr) {

  if (IsWeakDominated(pf, solPtr)) {
    return(pf)
  }

  wDom <- GetWeakDominated(pf, solPtr)
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
  return(pf)
}

HyperVolume <- function(pf) {
  hv <- dominated_hypervolume(pf$Getnorm(pf$ptrs),
                              matrix(1, pf$arch$dim, 1) * 1.1)
  return(hv)
}

IsWeakDominated <- function(pf, solPtr) {
  flag <- any(apply(FixRepmat(pf$arch$scores[solPtr, ], length(pf$ptrs), 1) >= pf$arch$scores[pf$ptrs, ], 1, all))
  return(flag)
}

GetWeakDominated <- function(pf, solPtr) {
  wDom <- apply(pf$arch$scores[pf$ptrs, ] >= FixRepmat(pf$arch$scores[solPtr, ], length(pf$ptrs), 1), 1, all)
  return(wDom)
}

UpdateMinMaxSc <- function(pf) {
  temp <- pf$arch$scores[pf$ptrs, ]
  if (is.matrix(temp)) {
    pf$scmax <- apply(temp, 2, max)
    pf$scmin <- apply(temp, 2, min)
  } else {
   pf$scmax <- temp
   pf$scmin <- temp
  }
  return(pf)
}

SetMinMaxSc <- function(pf, scmax, scmin) {
  pf$scmax <- scmax
  pf$scmin <- scmin
  return(pf)
}

GetNorm <- function(pf, solInd) {
  cnorm <- (pf$arch$scores[solInd, ] - FixRepmat(pf$scmin, length(solInd), 1)) /
    FixRepmat(pf$scmax - pf$scmin, length(solInd), 1)
  return(cnorm)
}
