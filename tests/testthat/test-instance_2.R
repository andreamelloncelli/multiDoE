options(digits = 10)

# setting
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
#### test MSOpt e Score ####

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

test_that("Score works",{expect_equal(Score(msopt, example),
                                      c(0.39519450275, 0.33394906171,
                                        0.08396341234, 0.11292875087,
                                        0.08464493702, 0.08167875103),
                                      tolerance = 0.0000000001
)
})

#### test MSSearch Single Crit ####
set.seed(13)
criteria <- "Id"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)

load("C:\\Users\\Francesca\\Desktop\\Rtesi\\multiDoE\\mssearch1_i2.RData")

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100),
               list("optsol" = mssearch1$optsol,
                    "optsc" = 0.3310033528,
                    "feval" = 271923,
                    "trend" = mssearch1$trend
               )
  )
})

#### test TPLSearch ####
set.seed(345)
criteria <-  c('I', 'Id', 'D')

lCrit <- length(criteria)
iters <- 10 * lCrit
restarts <- 100
restInit <- 2
i = 1

load("tpls3_i2.RData")
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





