options(digits = 10)
set.seed(123)

# setting
facts <- list(1:3, 4:6)
units <- list(12, 4)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "quadratic"

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
example <- matrix(c( 0,  1,  1, -1,  0,  0,
                     0,  1,  1,  0,  0,  1,
                     0,  1,  1,  0,  1, -1,
                     0,  1,  1,  1, -1,  0,
                     0, -1,  0,  0, -1,  0,
                     0, -1,  0,  0,  0,  0,
                     0, -1,  0,  1,  1,  1,
                     0, -1,  0, -1,  0, -1,
                     1, -1, -1,  1, -1,  1,
                     1, -1, -1,  1,  1, -1,
                     1, -1, -1, -1, -1, -1,
                     1, -1, -1, -1,  1,  1,
                     -1,  1, -1, -1,  1,  1,
                     -1,  1, -1,  1, -1,  1,
                     -1,  1, -1, -1, -1, -1,
                     -1,  1, -1,  1,  1, -1,
                     1,  0,  1,  0,  1,  0,
                     1,  0,  1, -1, -1,  1,
                     1,  0,  1,  0, -1, -1,
                     1,  0,  1,  1,  0,  0,
                     0,  0, -1,  0,  0,  0,
                     0,  0, -1,  1,  1,  0,
                     0,  0, -1, -1, -1,  1,
                     0,  0, -1,  0,  0, -1,
                     -1, -1, -1,  1,  1,  1,
                     -1, -1, -1,  1, -1, -1,
                     -1, -1, -1, -1,  1, -1,
                     -1, -1, -1, -1, -1,  1,
                     1,  1,  0,  1,  0, -1,
                     1,  1,  0,  1,  1,  1,
                     1,  1,  0,  0, -1,  0,
                     1,  1,  0, -1,  1, -1,
                     0,  0,  0,  1, -1, -1,
                     0,  0,  0,  0,  0,  1,
                     0,  0,  0, -1,  0,  0,
                     0,  0,  0,  0,  0,  0,
                     -1,  0,  0, -1,  1, -1,
                     -1,  0,  0,  0,  0,  0,
                     -1,  0,  0,  1,  0,  0,
                     -1,  0,  0,  0, -1,  0,
                     -1, -1,  1,  1,  1, -1,
                     -1, -1,  1,  1, -1,  1,
                     -1, -1,  1, -1,  1,  1,
                     -1, -1,  1, -1, -1, -1,
                     -1,  1,  1,  1,  1,  1,
                     -1,  1,  1,  0,  0, -1,
                     -1,  1,  1, -1,  1,  0,
                     -1,  1,  1, -1, -1,  1),
                  ncol = 6, byrow = T)

#### test ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1:3, 4:6),
                    "nfacts" = 6,
                    "nstrat" = 2,
                    "units" = list(12, 4),
                    "runs" = 48,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1),
                                   c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3, 3, 3, 3),
                    "Vinv" = t(solve(diag(48) + 1 * kronecker(diag(12), matrix(1, 4, 4)))),
                    "model"  = "quadratic",
                    "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
                    "ncrit" = 6,
                    "M" = Mq,
                    "M0" = M0q,
                    "W" = Wq)
               )
  }
)

test_that("Score works",
          {expect_equal(Score(msopt, example),
                        c(0.822843249300000, 0.719189332100000, 0.087740799100000,
                          0.176822808600000, 0.087184909200000, 0.119542430200000),
                        tolerance = 0.00001
          )
          }
)



