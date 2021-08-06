options(digits = 10)
set.seed(123)

# setting
facts <- list(1, 2:4, 5)
units <- list(21, 2, 1)
levels <- c(3, 4, 2, 2, 2)
etas <- list(1, 1)
criteria <- c('I','Id', 'D', 'A', 'Ds', 'As')
model <- "interaction"

# M
k <- length(unlist(facts))
k2 <- k * (k - 1) / 2
Mq <- rbind(
  cbind(1, t(integer(k)), t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
)

# M0
M0q <- rbind(
  cbind(1, t(integer(k)), t(integer(k2))),
  cbind(integer(k), diag(k) / 3, matrix(0, k, k2)),
  cbind(integer(k2), matrix(0, k2, k), diag(k2) / 9)
)
M0q[1, ] <- 0
M0q[, 1] <- 0

# W
nfacts <-length(unlist(facts))
w <- t(rep(1, nfacts + nfacts * (nfacts - 1) / 2))
a <- length(w / sum(w))
Wq <- c(w / sum(w)) * diag(a)

# Vinv
V <- diag(42)
for (i in 1:2) {
  if (i + 1 > length(units)) {
    ones_shape <- 1
  } else {
    ones_shape <- prod(unlist(units)[(i + 1):length(units)])
  }
  V <- V + etas[[1]] * kronecker(diag(prod(unlist(units)[1:i])),
                                matrix(1, ones_shape, ones_shape))
}
Vinv <- t(solve(V))

# msopt
msopt <- MSOpt(facts, units, levels, etas, criteria, model)

# example
example <- matrix(c( 1,  1,  1,  0,  0,
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

#### test ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1, 2:4, 5),
                    "nfacts" = 5,
                    "nstrat" = 3,
                    "units" = list(21, 2, 1),
                    "runs" = 42,
                    "etas" = list(1,1),
                    "avlev" = list(c(-1, 0 , 1), c(-1.0000000000, -0.3333333333,
                                                   0.3333333333, 1.0000000000),
                                   c(-1, 1), c(-1, 1), c(-1, 1)),
                    "levs" = c(3, 4, 2, 2, 2),
                    "Vinv" =  Vinv,
                    "model"  = "interaction",
                    "crit" =  c('I','Id', 'D', 'A', 'Ds', 'As'),
                    "ncrit" = 6
                    ,
                    "M" = Mq,
                    "M0" = M0q,
                    "W" = Wq))
}
)

test_that("Score works",
          {expect_equal(Score(msopt, example),
                        c( 0.413033963800000, 0.316643309000000, 0.114144889200000,
                           0.119475527300000, 0.115531263300000, 0.121014518800000),
                        tolerance = 0.00001)
          }
)




