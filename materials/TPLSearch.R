
# function TPLSearch ---------------------------

# library(bazar)
# library(scatterplot3d)
# library(vrmlgen)

TPLSearch <- function(facts, units, criteria, model, ...) {

  varargin <- list(...)

  tplsearch <- list()     # contiene ar, stats
  tplsearch$ar <- list()
  tplsearch$stats <- 0

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
  if (nargs() > 5){          # if aggiunto
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
  # lista di lunghezza 6 con 6 soluzioni iniziali [[i]] matrice
   initScores <- matrix(0, nCrit, nCrit) # matrice 6 x 6
   totFEval <- 0 # numero
  #


   o = 0

   for (i in 1:nCrit) {
     a <- t(c(numeric(i - 1), 1, numeric(nCrit - i)))   # alpha
     mssearch <- MSSearch(mso, a, "Restarts", restInit)

     initSol[[i]] <- mssearch$optsol
     initScores[i, ] <- mssearch$optsc
     totFEval <- totFEval + mssearch$feval
   }

   print("il primo vettore è lungo così")

############################################ x confronto MATLAB
    uno <- as.matrix(read.csv2("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\initSol1.txt",
                                sep = ",", header = F, colClasses = c("double", "double", "double", "double")), ncol = 4, byrow = T)
    due <- as.matrix(read.csv2("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\initSol2.txt",
                                sep = ",", header = F, colClasses = c("double", "double", "double", "double")), ncol = 4, byrow = T)
    tre <- as.matrix(read.csv2("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\initSol3.txt",
                                sep = ",", header = F, colClasses = c("double", "double", "double", "double")), ncol = 4, byrow = T)
    qua <- as.matrix(read.csv2("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\initSol4.txt",
                                sep = ",", header = F, colClasses = c("double", "double", "double", "double")), ncol = 4, byrow = T)
    cin <- as.matrix(read.csv2("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\initSol5.txt",
                                sep = ",", header = F, colClasses = c("double", "double", "double", "double")), ncol = 4, byrow = T)
    sei <- as.matrix(read.csv2("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\initSol6.txt",
                               sep = ",", header = F, colClasses = c("double", "double", "double", "double")), ncol = 4, byrow = T)
    initSol2 <- list(uno, due, tre, qua, cin, sei)
    initScores2 <- as.matrix(read.csv2("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\initScores.txt",
                            sep = ",",  dec = ".", header = F, colClasses = c("double", "double", "double", "double", "double", "double")))
    totFEval2 <- as.numeric(read.csv("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\totFEval.txt",
                          sep = ",", header = F))
    initSol <- initSol2
    totFEval <- totFEval2
    initScores <- initScores2
############################################


  ar <- Archive(nCrit, scalarizations + nCrit)

  for (i in 1:nCrit) {
    ar <- Add(ar, initSol[[i]], initScores[i, ])
  }

  pf <- PFront(ar)

  if (length(tplsearch) > 1) {
    stats <- totFEval
  }

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
    norms <- c(pf$scmin, pf$scmax - pf$scmin)
    # r <- matrix(runif(nCrit - 1), 1, nCrit - 1)
    r <- as.matrix(read.csv("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\r.txt", sep = ",", header = F))
    alpha <- cbind(r, 1 - sum(r))
    # alpha <- as.matrix(read.csv("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\alpha.txt", sep = ",", header = F))
    # start <- pf$arch$solutions[pf$ptrs[sample(1, 1:length(pf$ptrs))]]
    start <- read.csv("C:\\Users\\Francesca\\Desktop\\multiDoE_zip\\start.txt", sep = ",", header = F)
    start <- list(as.matrix(start))


    print("entro nel secondo MSSearch")
    print(i)
    newMssearch <- MSSearch(mso, alpha, "Start", start, "Normalize", norms)
    newSol <- newMssearch$optsol
    newScore <- newMssearch$optsc
    newFeval <- newMssearch$feval

    pf$arch <- Add(pf$arch, newSol, newScore)
    pf <- Add_PF(pf, pf$arch$nsols)

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
      stats <- c(stats, newFeval + stats[which.max(stats)])
    }
  }

  tplsearch$ar <- pf$arch
  tplsearch$stats <- stats
  return(tplsearch)
}
