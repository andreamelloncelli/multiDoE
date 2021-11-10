setwd("C:/Users/Francesca/Desktop/new1/multiDoE")
options(digits = 10)

# setting
facts <- list(1:2, 3:4)
units <- list(12, 4)
levels <- c(4, 4, 4, 2)
etas <- list(1, 1)
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

example <- matrix(c(  1, -1,  1,  0,
                      1, -1, -1,  1,
                      1, -1, -1,  1,
                      1, -1, -1,  1,
                      1, -1,  1,  0,
                      1, -1, -1,  1,
                      1, -1, -1,  1,
                      1, -1, -1,  1,
                      1,  1,  1,  1,
                      1,  1,  0,  1,
                      1,  1, -1,  0,
                      1,  1,  1, -1,
                      1,  1,  1,  1,
                      1,  1,  0,  1,
                      1,  1, -1,  0,
                      1,  1,  1, -1,
                      -1, -1,  0,  1,
                      -1, -1,  1, -1,
                      -1, -1,  1, -1,
                      -1, -1,  0,  1,
                      0, -1,  0, -1,
                      0, -1,  1,  1,
                      0, -1, -1, -1,
                      0, -1, -1,  0,
                      1,  1,  0,  1,
                      1,  1,  1, -1,
                      1,  1,  1,  1,
                      1,  1, -1,  0,
                      -1,  1, -1,  1,
                      -1,  1, -1, -1,
                      -1,  1,  1, -1,
                      -1,  1,  1,  1,
                      0, -1,  0, -1,
                      0, -1,  1,  1,
                      0, -1, -1,  0,
                      0, -1, -1, -1,
                      0, -1,  0, -1,
                      0, -1,  1,  1,
                      0, -1, -1, -1,
                      0, -1, -1,  0,
                      0,  0,  1,  0,
                      0,  0,  1,  0,
                      0,  0,  1,  0,
                      0,  0,  1,  0,
                      -1,  1,  1, -1,
                      -1,  1,  1,  1,
                      -1,  1, -1, -1,
                      -1,  1, -1,  1),
                  ncol = 4, byrow = T)
#### test MSOpt e Score ####

test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1:2, 3:4),
                    "nfacts" = 4,
                    "nstrat" = 2,
                    "units" = list(12, 4),
                    "runs" = 48,
                    "etas" = list(1, 1),
                    "avlev" = list(c(-1.0000000000, -0.3333333333, 0.3333333333, 1.0000000000),
                                   c(-1.0000000000, -0.3333333333, 0.3333333333, 1.0000000000),
                                   c(-1.0000000000, -0.3333333333, 0.3333333333, 1.0000000000),
                                   c(-1, 1)),
                    "levs" = c(4, 4, 4, 2),
                    "Vinv" = t(solve(diag(48) + 1 * kronecker(diag(12), matrix(1, 4, 4)))),
                    "model"  = 'quadratic',
                    "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
                    "ncrit" = 6,
                    "M" = Mq,
                    "M0" = M0q,
                    "W" = Wq
               )
  )
})

test_that("Score works",{expect_equal(Score(msopt, example),
                                      c(1.205934492791982, 0.734833087126215,
                                        0.108876783413959, 0.406282668029913,
                                        0.109221258351505, 0.168837086116107),
                                      tolerance = 0.0000000001)}
          )

#### test MSSearch Single Crit ####
set.seed(13)
criteria <- "As"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100),
               list("optsol" = matrix(0, 48, 4),
                    "optsc" = as.matrix(Inf),
                    "feval" = 26500,
                    "trend" = rep(Inf, 100)
               )
  )
})



#### test TPLSearch ####
set.seed(3455)
criteria <-  c('D', 'Aa')
options(digits = 10)

facts <- list(1:2, 3:4)
units <- list(12, 4)
levels <- c(4, 4, 4, 2)
etas <- list(1, 1)
criteria <- c('D', 'A')
model <- "main"

lCrit <- length(criteria)
iters <- 10 * lCrit
restarts <- 100
restInit <- 2
i = 70

file_name <- here::here("tests/testthat/test_data/tpls_i5.Rds")
tpls <- readRDS(file = file_name)
ar <- tpls$ar
stats <- tpls$stats
megaAR <- tpls$megaAR

test_that("runTPLSearch works", {
  expect_equal(runTPLS(facts,units, criteria, model, iters, "Etas", etas,
                       "Levels", levels, "Restarts", restarts, "RestInit",
                       restInit, "RngSeed", i),
               list("ar" = ar, "stats" = stats, "megaAR" = megaAR)
  )
}
)
