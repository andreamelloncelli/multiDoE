setwd(here::here())

options(digits = 15)

####################################################################
# setting interaction ####
facts <- list(1, 2:5)
units <- list(21, 2)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "interaction"

# M ####
k <- length(unlist(facts))
k2 <- k * (k - 1) / 2
Mi <- rbind(
  cbind(1, t(integer(k)), t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
)

# M0 ####
M0i <- rbind(
  cbind(1, t(integer(k)), t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
)
M0i[1, ] <- 0
M0i[, 1] <- 0

# W ####
nfacts <- length(unlist(facts))
w <- t(rep(1, nfacts + nfacts * (nfacts - 1) / 2))
a <- length(w / sum(w))
Wi <- c(w / sum(w)) * diag(a)

# msopt ####
msopt <- MSOpt(facts, units, levels, etas, criteria, model)

### test MSOpt e Score interaction: OK ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1, 2:5),
                    "nfacts" = 5,
                    "nstrat" = 2,
                    "units" = list(21, 2),
                    "runs" = 42,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1),
                                   c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3, 3, 3),
                    "Vinv" = t(solve(diag(42) + 1 * kronecker(diag(21), matrix(1, 2, 2)))),
                    "model"  = 'interaction',
                    "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
                    "ncrit" = 6,
                    "M" = Mi,
                    "M0" = M0i,
                    "W" = Wi
               )
  )
})

test_that("Score works", {
  # TODO the function Score does not work properly
  expect_equal(Score(msopt, example),
               c(0.257854563945983, 0.185759114983032,
                 0.065379595918449, 0.069753041857244,
                 0.064995045171766, 0.069596881383530),
               tolerance = 0.0000000001
  )
})




####################################################################
# setting main ####
facts <- list(1, 2:5)
units <- list(21, 2)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "main"

# M ####
k <- length(unlist(facts))
k2 <- k * (k - 1) / 2
Mm <- rbind(
  cbind(1, t(integer(k))),
  cbind(integer(k), diag(k) / 3)
)

# M0 ####
M0m <- rbind(
  cbind(1, t(integer(k))),
  cbind(integer(k), diag(k) / 3)
)
M0m[1, ] <- 0
M0m[, 1] <- 0

# W ####
nfacts <- length(unlist(facts))
w <- t(rep(1, msopt$nfacts))
a <- length(w / sum(w))
Wm <- c(w / sum(w)) * diag(a)

# msopt ####
msopt <- MSOpt(facts, units, levels, etas, criteria, model)

### test MSOpt e Score main: OK ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1, 2:5),
                    "nfacts" = 5,
                    "nstrat" = 2,
                    "units" = list(21, 2),
                    "runs" = 42,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1),
                                   c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3, 3, 3),
                    "Vinv" = t(solve(diag(42) + 1 * kronecker(diag(21), matrix(1, 2, 2)))),
                    "model"  = 'main',
                    "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
                    "ncrit" = 6,
                    "M" = Mm,
                    "M0" = M0m,
                    "W" = Wm
               )
  )
})

test_that("Score works",{
  # TODO fix Score
  expect_equal(Score(msopt, example),
               c(0.173111218905576, 0.101633636740146,
                 0.057209958841188, 0.062729748730978,
                 0.054725730766041, 0.060980182044087),
               tolerance = 0.0000000001
  )
})






####################################################################
# setting quadratic ####
facts <- list(1, 2:5)
units <- list(21, 2)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "quadratic"

# M ####
k <- length(unlist(facts))
k2 <- k * (k - 1) / 2
Mq <- rbind(
  cbind(1, t(integer(k)), t(rep(1, k)) / 3, t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k), matrix(0, k, k2)),
  cbind(rep(1, k) / 3, matrix(0, k, k),
        (4 * diag(k) + 5 * matrix(1, k, k)) / 45, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), matrix(0, k2, k), diag(k2) / 9)
)

# M0 ####
M0q <- rbind(
  cbind(1, t(integer(k)), t(rep(1, k)) / 3, t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k), matrix(0, k, k2)),
  cbind(rep(1, k) / 3, matrix(0, k, k),
        (4 * diag(k) + 5 * matrix(1, k, k)) / 45, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), matrix(0, k2, k), diag(k2) / 9)
)
M0q[1, ] <- 0
M0q[, 1] <- 0

# W ####
nfacts <- length(unlist(facts))
w <- c(
  rep(1, nfacts),
  rep(1, nfacts) / 4,
  rep(1, nfacts * (nfacts - 1) / 2)
)
a <- length(w / sum(w))
Wq <- c(w / sum(w)) * diag(a)

# msopt ####
msopt <- MSOpt(facts, units, levels, etas, criteria, model)

# example ####
example <- matrix(c(1,  1,  1,  0,  0,
                    1, -1, -1,  1, -1,
                    1,  0, -1,  0,  0,
                    1,  1,  1, -1,  1,
                    1, -1,  1,  1, -1,
                    1, -1, -1, -1,  1,
                    0,  1, -1,  0,  0,
                    0,  1,  1,  1, -1,
                    0, -1, -1,  0,  0,
                    0,  0,  0,  1,  1,
                    0,  0,  1, -1,  0,
                    0, -1,  0,  0,  1,
                    0,  0,  1,  0,  1,
                    0, -1,  0, -1,  0,
                    0,  1, -1, -1,  1,
                    0,  0,  0,  0,  0,
                    -1,  0,  0,  1,  0,
                    -1,  1,  1, -1, -1,
                    1, -1,  1,  1,  1,
                    1,  1, -1,  1, -1,
                    0,  0,  0,  0,  0,
                    0, -1,  1, -1,  1,
                    -1, -1,  1, -1,  1,
                    -1,  1, -1,  1, -1,
                    0,  1, -1, -1, -1,
                    0,  0,  0,  0,  0,
                    -1,  0,  0,  0, -1,
                    -1,  1,  1,  1,  1,
                    0, -1, -1,  1,  0,
                    0,  0,  0,  0, -1,
                    -1, -1,  1,  1, -1,
                    -1,  0, -1, -1,  1,
                    0,  0, -1, -1, -1,
                    0,  1,  0,  0,  1,
                    -1, -1, -1,  1,  1,
                    -1,  1,  0, -1,  0,
                    1, -1,  1, -1, -1,
                    1,  0,  0,  1,  0,
                    1,  1,  0, -1, -1,
                    1,  1, -1,  1,  1,
                    -1,  0,  1,  0,  0,
                    -1, -1, -1, -1, -1),
                  ncol = 5, byrow = T)
#### test MSOpt e Score quadratic: OK ####

test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1, 2:5),
                    "nfacts" = 5,
                    "nstrat" = 2,
                    "units" = list(21, 2),
                    "runs" = 42,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1),
                                   c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3, 3, 3),
                    "Vinv" = t(solve(diag(42) + 1 * kronecker(diag(21), matrix(1, 2, 2)))),
                    "model"  = 'quadratic',
                    "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
                    "ncrit" = 6,
                    "M" = Mq,
                    "M0" = M0q,
                    "W" = Wq
               )
  )
})

test_that("Score works",{
  expect_equal(Score(msopt, example),
               c(0.395194502748104, 0.333949061712227,
                 0.083963412344681, 0.112928750871800,
                 0.084644937016109, 0.081678751029379),
               tolerance = 0.0000000001
  )
})



####################################################################
#### test MSSearch Single Crit ("Id") + Restarts: OK ####
set.seed(13)
criteria <- "Id"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)
file_name <- here::here("tests/testthat/test_data/mss1_i2.Rds")
mssearch1 <- readRDS(file = file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100),
               list("optsol" = mssearch1$optsol,
                    "optsc" = 0.331003352814033,
                    "feval" = 271923,
                    "trend" = mssearch1$trend
               )
  )
})

#### test MSSearch Single Crit ("Id") + Restarts + Start: OK ####
set.seed(13)
criteria <- "Id"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)
file_name <- here::here("tests/testthat/test_data/mss1sol_i2.Rds")
mssearch1 <- readRDS(file = file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100, "Start", example),
               list("optsol" = mssearch1$optsol,
                    "optsc" = 0.331977510344595,
                    "feval" =  38659,
                    "trend" = mssearch1$trend
               )
  )
})

#### test TPLSearch ####
# set.seed(345)
# criteria <-  c('I', 'Id', 'D')
#
# lCrit <- length(criteria)
# iters <- 10 * lCrit
# restarts <- 100
# restInit <- 2
# i = 1
#
# load("tpls3_i2.RData")
# ar <- tpls$ar
# stats <- tpls$stats
# megaAR <- tpls$megaAR
#
# test_that("runTPLSearch works", {
#   expect_equal(runTPLS(facts,units, criteria, model, iters, "Etas", etas,
#                        "Levels", levels, "Restarts", restarts, "RestInit",
#                        restInit, "RngSeed", i),
#                list("ar" = ar, "stats" = stats, "megaAR" = megaAR)
#   )
# }
# )





