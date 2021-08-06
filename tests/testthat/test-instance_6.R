options(digits = 10)
library(multiDoE)

# setting
facts <- list(c(), 1:3)
units <- list(7, 4)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- 'quadratic'

# M
k <- length(unlist(facts))
k2 <- k * (k - 1) / 2
Mq <- rbind(
  cbind(1, t(integer(k)), t(rep(1, k)) / 3, t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k), matrix(0, k, k2)),
  cbind(rep(1, k) / 3, matrix(0, k, k),
        (4 * diag(k) + 5 * matrix(1, k, k)) / 45, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), matrix(0, k2, k), diag(k2) / 9)
)

# M0
M0q <- rbind(
  cbind(1, t(integer(k)), t(rep(1, k)) / 3, t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k), matrix(0, k, k2)),
  cbind(rep(1, k) / 3, matrix(0, k, k),
        (4 * diag(k) + 5 * matrix(1, k, k)) / 45, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), matrix(0, k2, k), diag(k2) / 9)
)
M0q[1, ] <- 0
M0q[, 1] <- 0

# W
nfacts <-length(unlist(facts))
w <- c(
  rep(1, nfacts),
  rep(1, nfacts) / 4,
  rep(1, nfacts * (nfacts - 1) / 2)
)
a <- length(w / sum(w))
Wq <- c(w / sum(w)) * diag(a)

# msopt
msopt <- MSOpt(facts, units, levels, etas, criteria, model)

# example
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

#### test ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(c(), 1:3),
                    "nfacts" = 3,
                    "nstrat" = 2,
                    "units" = list(7, 4),
                    "runs" = 28,
                    "etas" = list(1),
                    "avlev" = list(c(-1.0000000000, 0, 1.0000000000),
                                   c(-1.0000000000, 0, 1.0000000000),
                                   c(-1.0000000000, 0, 1.0000000000)),
                    "levs" = c(3, 3, 3),
                    "Vinv" = t(solve(diag(28) + 1 * kronecker(diag(7), matrix(1, 4, 4)))),
                    "model"  = "quadratic",
                    "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
                    "ncrit" = 6,
                    "M" = Mq,
                    "M0" = M0q,
                    "W" = Wq
               )
  )
}
)

test_that("Score works",
          {expect_equal(Score(msopt, example),
                        c(0.346163804300000, 0.177527693800000, 0.10306494420000,
                          0.128787948000000, 0.096959114100000, 0.086644705600000),
                        tolerance = 0.00001)
          }
)



