library(multiDoE)
options(digits = 10)

# setting
facts <- list(1:2, 3:4)
units <- list(12, 4)
levels <- c(4, 4, 4, 2)
etas <- list(1, 1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "main"

# M
k <- length(unlist(facts))
k2 <- k * (k - 1) / 2
Mq <- rbind(
  cbind(1, t(integer(k))),
  cbind(integer(k), diag(k) / 3)
)

# M0
M0q <- rbind(
  cbind(1, t(integer(k))),
  cbind(integer(k), diag(k) / 3)
)
M0q[1, ] <- 0
M0q[, 1] <- 0

# W
nfacts <-length(unlist(facts))
w <- t(rep(1, nfacts))
a <- length(w / sum(w))
Wq <- c(w / sum(w)) * diag(a)

# msopt
msopt <- MSOpt(facts, units, levels, etas, criteria, model)

# example
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

#### test ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(1:2, 3:4),
                    "nfacts" = 4,
                    "nstrat" = 2,
                    "units" = list(12, 4),
                    "runs" = 48,
                    "etas" = list(1,1),
                    "avlev" = list(c(-1.0000000000, -0.3333333333, 0.3333333333, 1.0000000000),
                                   c(-1.0000000000, -0.3333333333, 0.3333333333, 1.0000000000),
                                   c(-1.0000000000, -0.3333333333, 0.3333333333, 1.0000000000),
                                   c(-1, 1)),
                    "levs" = c(4, 4, 4, 2),
                    "Vinv" = t(solve(diag(48) + 1 * kronecker(diag(12), matrix(1, 4, 4)))),
                    "model"  = "main",
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
                        c(0.224788195600000, 0.114820240900000, 0.072016241900000,
                          0.090885735500000, 0.065668315500000, 0.086115180700000),
                          tolerance = 0.00001)
          }
)


