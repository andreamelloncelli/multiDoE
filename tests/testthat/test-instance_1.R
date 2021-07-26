options(digits = 10)
set.seed(123)

facts <- list(1, 2:5)
units <- list(6, 5)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "quadratic"

# M
k <- 5
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
w <- c(
  rep(1, 5),
  rep(1, 5) / 4,
  rep(1, 5 * (5 - 1) / 2)
)
a <- length(w / sum(w))
Wq <- c(w / sum(w)) * diag(a)


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
          {expect_equal(Score(list("facts" = list(1, 2:5),
                                   "nfacts" = 5,
                                   "nstrat" = 2,
                                   "units" = list(6, 5),
                                   "runs" = 30,
                                   "etas" = list(1),
                                   "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1),
                                                  c(-1, 0, 1), c(-1, 0, 1),
                                                  c(-1, 0, 1)),
                                   "levs" = c(3, 3, 3, 3, 3),
                                   "Vinv" = t(solve(diag(30) + 1 * kronecker(diag(6), matrix(1, 5, 5)))),
                                   "model"  = "quadratic",
                                   "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
                                   "ncrit" = 6,
                                   "M" = Mq,
                                   "M0" = M0q,
                                   "W" = Wq
                                   ),
                              settings = example
                              ),
  c(0.7476857097, 0.5437606070, 0.1016563124, 0.1826613447, 0.0982742118, 0.1046020401),
  tolerance = 0.00001
  )
}
)







