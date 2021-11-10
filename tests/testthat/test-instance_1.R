setwd("C:/Users/Francesca/Desktop/new1/multiDoE")
options(digits = 10)

# setting
facts <- list(1, 2:5)
units <- list(6, 5)
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
example <- matrix(c( 0,  1,  0,  0,  0,
                     0,  0,  0,  1,  1,
                     0,  0, -1, -1,  0,
                     0,  0,  1,  0, -1,
                     0, -1,  0,  0,  0,
                     0,  0, -1,  0,  0,
                     0,  0,  1,  0,  0,
                     0, -1,  0,  0,  0,
                     0,  0,  0,  1, -1,
                     0,  1,  0, -1, -1,
                     -1, -1,  1, -1,  1,
                     -1,  1, -1, -1,  1,
                     -1, -1, -1,  1,  1,
                     -1,  1,  1,  1,  1,
                     -1,  0,  0,  0, -1,
                     1, -1, -1,  1, -1,
                     1,  0,  0, -1,  0,
                     1,  1,  1,  1, -1,
                     1, -1,  1,  1,  1,
                     1,  1, -1,  1,  1,
                     -1,  0,  0,  0,  1,
                     -1,  1, -1,  1, -1,
                     -1,  1,  1, -1,  0,
                     -1, -1, -1, -1, -1,
                     -1, -1,  1,  1, -1,
                     1,  1,  1, -1,  1,
                     1, -1,  1, -1, -1,
                     1,  1, -1,  0, -1,
                     1, -1, -1, -1,  1,
                     1,  0,  0,  1,  0),
                  ncol = 5, byrow = T)

#### test MSOpt e Score: OK ####

test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1, 2:5),
                    "nfacts" = 5,
                    "nstrat" = 2,
                    "units" = list(6, 5),
                    "runs" = 30,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1),
                                   c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3, 3, 3),
                    "Vinv" = t(solve(diag(30) + 1 * kronecker(diag(6), matrix(1, 5, 5)))),
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
                        c(0.747685709663110, 0.543760607038438, 0.101656312421005,
                          0.182661344673558, 0.098274211834103, 0.104602040115434),
                        tolerance = 0.0000000001
                        )
            })
#
#### test MSSearch Single Crit ("I") + Restarts: OK ####
set.seed(1)
criteria <- "I"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)
file_name <- here::here("tests/testthat/test_data/mss1_i1.RDs")
mssearch1 <- readRDS(file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100),
               list("optsol" = mssearch1$optsol,
                    "optsc" = 0.7602652588,
                    "feval" = 176361,
                    "trend" = mssearch1$trend
                    )
               )
  })

#### test MSSearch Single Crit ("I") + Restarts + Start: oK ####
file_name <- here::here("tests/testthat/test_data/mss1sol_i1.RDs")
mssearch1 <- readRDS(file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100, "Start", example),
               list("optsol" = mssearch1$optsol,
                    "optsc" = mssearch1$optsc,
                    "feval" = mssearch1$feval,
                    "trend" = mssearch1$trend
               )
  )
})

#### test MSSearch Single Crit ("I") + Restarts + Start + Normalize: OK ####
file_name <- here::here("tests/testthat/test_data/mss1norm_i1.Rds")
mssearch1 <- readRDS(file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100, "Start", example,
                        "Normalize", c(0.1, 2)),
               list("optsol" = mssearch1$optsol,
                    "optsc" = mssearch1$optsc,
                    "feval" = mssearch1$feval,
                    "trend" = mssearch1$trend
               )
  )
})




#### test MSSearch Multi Crit (tutti) + Restarts + Start: OK ####
file_name <- here::here("tests/testthat/test_data/mssMsol_i1.Rds")
mssearch1 <- readRDS(file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt, rep(1/6, 6), "Restarts", 100, "Start", example),
               list("optsol" = mssearch1$optsol,
                    "optsc" = mssearch1$optsc,
                    "feval" = mssearch1$feval,
                    "trend" = mssearch1$trend
               )
  )
})

#### test TPLSearch ####
set.seed(36)

facts <- list(1, 2:5)
units <- list(6, 5)
levels <- 3
etas <- list(1)
criteria <-  c('I', 'Id', 'D')
model <- "quadratic"

lCrit <- length(criteria)
iters <- 10 * lCrit
restarts <- 100
restInit <- 2
i = 70

tpls <- readRDS(here::here("tests/testthat/test_data/tpls_i1.Rds"))

test_that("runTPLSearch works", {
  expect_equal(runTPLS(facts, units, criteria, model, iters,
                       "Etas", etas,
                       "Levels", levels,
                       "Restarts", restarts,
                       "RestInit", restInit,
                       "RngSeed", i),
               list("ar" = tpls$ar, "stats" = tpls$stats, "megaAR" = tpls$megaAR)
  )
  }
)





