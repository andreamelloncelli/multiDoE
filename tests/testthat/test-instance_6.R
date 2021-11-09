options(digits = 10)

# setting
facts <- list(c(), 1:3)
units <- list(7, 4)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- 'quadratic'

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
example <- matrix(c( 0,  0,  0,
                     1, -1, -1,
                     -1,  1, -1,
                     -1, -1,  1,
                     1, -1,  1,
                     0,  1,  0,
                     -1, -1, -1,
                     0,  0,  0,
                     0,  0,  0,
                     1,  0,  1,
                     0, -1, -1,
                     -1,  1,  0,
                     1,  1,  0,
                     0, -1,  1,
                     -1,  1,  1,
                     -1,  0, -1,
                     0, -1,  0,
                     0,  0,  1,
                     1,  1, -1,
                     -1,  0,  0,
                     -1, -1,  0,
                     0,  1,  1,
                     1,  0, -1,
                     0,  0,  0,
                     1,  1,  1,
                     0,  1, -1,
                     -1,  0,  1,
                     1, -1,  0),
                  ncol = 3, byrow = T)


#### test MSOpt e Score ####

test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(c(), 1:3),
                    "nfacts" = 3,
                    "nstrat" = 2,
                    "units" = list(7, 4),
                    "runs" = 28,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3),
                    "Vinv" = t(solve(diag(28) + 1 * kronecker(diag(7), matrix(1, 4, 4)))),
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
                                      c(0.346163804255849,
                                        0.177527693785788,
                                        0.103064944168469,
                                        0.128787947994373,
                                        0.096959114137211,
                                        0.086644705619155),
                                      tolerance = 0.0000000001)
}
)


#### test MSSearch Single Crit ("Ds") + Restarts ####
set.seed(13)
criteria <- "Ds"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)

# file_name <- here::here("tests/testthat/test_data/mss1_i6.Rds")
file_name <- here::here(file.path("tests", "testthat","test_data", "mss1_i6.Rds"))
mssearch1 <- readRDS(file = file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100),
               list("optsol" = mssearch1$optsol,
                    "optsc" = 0.08567553263,
                    "feval" =  87413,
                    "trend" = mssearch1$trend
               )
  )
})




#### test MSSearch Single Crit ("Ds") + Restarts + Start ####
set.seed(13)

criteria <- "Ds"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)
# file_name <- here::here("tests/testthat/test_data/mss1sol_i6.Rds")
file_name <- here::here(file.path("tests", "testthat","test_data", "mss1sol_i6.Rds"))
mssearch1 <- readRDS(file = file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100, "Start", example),
               list("optsol" = mssearch1$optsol,
                    "optsc" = 0.08780585259,
                    "feval" =  17245,
                    "trend" = rep(0.08780585259, 100)
               )
  )
})





#### test MSSearch Multi Crit ####
set.seed(13)
facts <- list(c(), 1:3)
units <- list(7, 4)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "quadratic"

msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)
file_name <- here::here("tests/testthat/test_data/mssM_i6.Rds")
mssearchm <- readRDS(file = file_name)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, rep(1/6, 6), "Restarts", 100, "Start", example,
                        "Normalize", c(rep(0.2,6), rep(0.5,6))),
               list("optsol" = mssearchm$optsol,
                    "optsc" = c(0.346163804255849, 0.177527693785788,
                                0.103064944168469, 0.128787947994373,
                                0.096959114137211, 0.086644705619155),
                    "feval" =  16900,
                    "trend" = mssearchm$trend
               ), tolerance = 0.0000000001
  )
})






set.seed(345)
options(digits = 10)

facts <- list(c(), 1:3)
units <- list(7, 4)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id')
model <- 'main'

lCrit <- length(criteria)
iters <- 10 * lCrit
restarts <- 100
restInit <- 2
i = 70

#file_name <- here::here("tests/testthat/test_data/tpls_i6.Rds")
file_name <- here::here(file.path("tests", "testthat","test_data", "tpls_i6.Rds"))
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
