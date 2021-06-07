
# function TPLSearch ---------------------------

# library(bazar)
# library(scatterplot3d)
# library(vrmlgen)

TPLSearch <- function(facts, units, criteria, model, ...) {

  varargin <- list(...)
  tplsearch <- list()     # contiene ar, stats

  nStrat <- length(facts)
  nCrit <- length(criteria)

  # default parameters
  restarts <- 1000
  levels <- 3
  etas <- matrix(1, 1, nStrat - 1)
  rngSeed <- 0
  startMode <- 1
  alphaMode <- 0
  restInit <- 50
  interact <- 0   # interactive è nome di funzione

  # optional parameters
  if (nargs() > 5){          # ho dovuto aggiungere l'if
    for (i in seq(1, nargs() - 5, 2)) {
      switch(varargin[[i]],
             "Restarts" = {
               restarts = varargin[[i+1]]
             },
             "Levels" = {
               levels = varargin[[i+1]]
             },
             "Etas" = {
               etas = varargin[[i+1]]
             },
             "RngSeed" = {
               rngSeed = varargin[[i+1]]
             },
             "StartMode" = {
               startMode = varargin[[i+1]]
             },
             "AlphaMode" = {
               alphaMode = varargin[[i+1]]
             },
             "RestInit" = {
               restInit = varargin[[i+1]]
             },
             "Interactive" = {
               interact = varargin[[i+1]]
             }
      )
    }
  }

  scalarizations <- restarts - (restInit * nCrit)
  mso <- MSOpt(facts, units, levels, etas, criteria, model)
  set.seed(rngSeed)

  # good quality solutions for the single objectives
  initSol <- vector(mode = "list", length = nCrit)
  initScores <- matrix(0, nCrit, nCrit)
  totFEval <- 0

  print("fin qui ok")


  for (i in 1:nCrit) {
    print('number of crit')
    print(i)
    a <- t(c(numeric(i - 1), 1, numeric(nCrit - i)))   # alpha
    mssearch <- MSSearch(mso, a, "Restarts", restInit)   # modifica con a[i]
    print(mso)
    initSol[[i]] <- mssearch$optsol
    initScores[i, ] <- mssearch$optsc
    totFEval <- totFEval + mssearch$feval
    print("optsc")
    print(mssearch$optsc)

  }
  print("before ar")
  ar <- Archive(nCrit, scalarizations + nCrit)

  print("before add for")
  for (i in 1:nCrit) {
    Add(ar, initSol[[i]], initScores[i, ])
  }
  print("before pf")

  pf <- PFront(ar)

  if (length(tplsearch) > 1) {
    stats <- totFEval
  }

  print("before interact")
  if (interact) {
    if (nCrit==2){
      plot(ar$scores[1:ar$nsols, 1], ar$scores[1:ar$nsols, 2], col = "blue")
      points(ar$scores[pf$ptrs, 1], ar$scores[pf$ptrs, 2], col = "red")
      #pause(Inf)
    }
    else if (nCrit==3){
      scatterplot3d(ar$score[1:ar$nsols, 1], ar$score[1:ar$nsols, 2],ar$score[1:ar$nsols, 3], color = "blue")
      point3d(ar$score[pf$ptrs, 1], ar$score[pf$ptrs, 2],ar$score[pf$ptrs, 3], color = "red")
      #pause(Inf)
    }
  }

  for (i in 1:scalarizations) {
    norms <- cbind(pf$scmin, pf$scmax - pf$scmin)
    r <- matrix(runif(nCrit - 1), 1, nCrit - 1)
    alpha <- cbind(r, 1 - sum(r))
    start <- ar$solutions[pf$ptrs(sample(1, 1:length(pf$ptrs)))]

    newMssearch <- MSSearch(mso, alpha, "Start", start, "Normalize", norms)
    newSol <- newMssearch$optsol
    newScore <- newMssearch$optsc
    newFeval <- newMssearch$feval

    Add(ar, newSol, newScore)
    Add_PF(pf, ar$nsols)

    if (interact) {
      if (nCrit==2){
        plot(ar$scores[1:ar$nsols, 1], ar$scores[1:ar$nsols, 2], col = "blue")
        points(ar$scores[pf$ptrs, 1], ar$scores[pf$ptrs, 2], col = "red")
        #pause(Inf)
      }
      else if (nCrit==3){
        scatterplot3d(ar$score[1:ar$nsols, 1], ar$score[1:ar$nsols, 2],ar$score[1:ar$nsols, 3], color = "blue")
        point3d(ar$score[pf$ptrs, 1], ar$score[pf$ptrs, 2],ar$score[pf$ptrs, 3], color = "red")
        #pause(Inf)
      }
    }
    if (length(tplsearch) > 1) {
      stats <- cbind(stats, newFeval + as.matrix(stats)[dim(stats)[1], 1])
    }
  }

  tplsearch$ar <- ar
  tplsearch$stats <- stats
  return(tplsearch)

}
