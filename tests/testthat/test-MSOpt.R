M0q <- rbind(
  cbind(1, t(integer(4)), t(rep(1, 4)) / 3, t(integer(6))),
  cbind(integer(4), diag(4) / 3, matrix(0, 4, 4), matrix(0, 4, 6)),
  cbind(rep(1, 4) / 3, matrix(0, 4, 4),
        (4 * diag(4) + 5 * matrix(1, 4, 4)) / 45, matrix(0, 4, 6)),
  cbind(integer(6), matrix(0, 6, 4), matrix(0, 6, 4), diag(6) / 9)
  )
M0q[1, ] <- 0
M0q[, 1] <- 0


Mq <- rbind(
  cbind(1, t(integer(4)), t(rep(1, 4)) / 3, t(integer(6))),
  cbind(integer(4), diag(4) / 3, matrix(0, 4, 4), matrix(0, 4, 6)),
  cbind(rep(1, 4) / 3, matrix(0, 4, 4), (4 * diag(4) + 5 * matrix(1, 4, 4)) / 45,
        matrix(0, 4, 6)),
  cbind(integer(6), matrix(0, 6, 4), matrix(0, 6, 4), diag(6) / 9)
)


w <-  c(
  rep(1, 4),
  rep(1, 4) / 4,
  rep(1, 4 * (4 - 1) / 2)
  )
a <- length(w / sum(w))
Wq <- c(w / sum(w)) * diag(a)


test_that("MSOpt works", {
  expect_equal(MSOpt(list(1:2, 3:4), list(12, 4), 3, list(1, 1),
                     c('I', 'Id', 'D', 'A', 'Ds', 'As'), 'quadratic'
                     ),
               list("facts" = list(1:2, 3:4),
                    "nfacts" = 4,
                    "nstrat" = 2,
                    "units" = list(12, 4),
                    "runs" = 48,
                    "etas" = list(1, 1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3, 3),
                    "Vinv" = t(solve(diag(48) + 1 * kronecker(diag(12), matrix(1, 4, 4)))),
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


test_that("Score works", {
  expect_equal(Score(list(
    "facts" = list(1:2, 3:4),
    "nfacts" = 4,
    "nstrat" = 2,
    "units" = list(12, 4),
    "runs" = 48,
    "etas" = list(1, 1),
    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)),
    "levs" = c(3, 3, 3, 3),
    "Vinv" = t(solve(diag(48) + 1 * kronecker(diag(12), matrix(1, 4, 4)))),
    "model"  = "quadratic",
    "crit" = c('I', 'Id', 'D', 'A', 'Ds', 'As'),
    "ncrit" = 6,
    "M" = Mq,
    "M0" = M0q,
    "W" = Wq
    ),
    settings = example
    ),
    c(1.2059, 0.7348, 0.1089, 0.4063, 0.1092, 0.1688),
    tolerance = 0.0001
  )
}
)





