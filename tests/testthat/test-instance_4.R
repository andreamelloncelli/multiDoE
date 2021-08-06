options(digits = 10)
set.seed(123)

# setting
facts <- list(1:2, 3:4)
units <- list(10, 5)
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
example <- matrix(c( -1,  1, -1, -1,
                     -1,  1, -1,  1,
                     -1,  1,  1,  1,
                     -1,  1,  0,  0,
                     -1,  1,  1, -1,
                     0,  0, -1,  0,
                     0,  0,  1,  0,
                     0,  0,  0,  1,
                     0,  0,  0, -1,
                     0,  0,  0,  0,
                     0,  0,  0,  0,
                     0,  0,  0, -1,
                     0,  0,  0,  1,
                     0,  0, -1,  0,
                     0,  0,  1,  0,
                     1,  0,  1,  0,
                     1,  0, -1,  1,
                     1,  0,  0, -1,
                     1,  0, -1, -1,
                     1,  0,  0,  0,
                     -1,  0,  0, -1,
                     -1,  0, -1, -1,
                     -1,  0, -1,  1,
                     -1,  0,  0,  0,
                     -1,  0,  1,  0,
                     0,  1,  1, -1,
                     0,  1,  1,  1,
                     0,  1,  0,  0,
                     0,  1,  0,  1,
                     0,  1, -1,  0,
                     1, -1, -1, -1,
                     1, -1,  1,  1,
                     1, -1, -1,  1,
                     1, -1,  1, -1,
                     1, -1,  0,  0,
                     1,  1,  1,  1,
                     1,  1,  1, -1,
                     1,  1, -1,  1,
                     1,  1, -1, -1,
                     1,  1,  0,  0,
                     -1, -1,  1,  1,
                     -1, -1,  0,  0,
                     -1, -1, -1, -1,
                     -1, -1, -1,  1,
                     -1, -1,  1, -1,
                     0, -1,  1, -1,
                     0, -1,  0,  0,
                     0, -1, -1,  0,
                     0, -1,  0,  1,
                     0, -1,  1,  1),
                  ncol = 4 , byrow = T)

#### test ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1:2, 3:4),
                    "nfacts" = 4,
                    "nstrat" = 2,
                    "units" = list(10, 5),
                    "runs" = 50,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1),
                                   c(-1, 0, 1)),
                    "levs" = c(3, 3, 3, 3),
                    "Vinv" = t(solve(diag(50) + 1 * kronecker(diag(10), matrix(1, 5, 5)))),
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
                        c(0.511761111700000, 0.426046826000000, 0.097706485000000,
                          0.179635486900000, 0.096282604200000, 0.118533289500000),
                        tolerance = 0.00001
          )
          }
)




