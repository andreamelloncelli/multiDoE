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
                                      c(0.34616380426, 0.17752769379,
                                        0.10306494417, 0.12878794799,
                                        0.09695911414, 0.08664470562),
                                      tolerance = 0.0000000001)
}
)


#### test MSSearch Single Crit ####
set.seed(13)
criteria <- "Ds"
msopt1 <- MSOpt(facts, units, levels, etas, criteria, model)

load("C:\\Users\\Francesca\\Desktop\\Rtesi\\multiDoE\\mssearch1_i6.RData")

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt1, 1, "Restarts", 100),
               list("optsol" = mssearch1$optsol,
                    "optsc" = 0.08567553263,
                    "feval" =  87413,
                    "trend" = mssearch1$trend
               )
  )
})



