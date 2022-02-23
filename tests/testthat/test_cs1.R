# instance 2 paper 2016
setwd(here::here())
options(digits = 10)

# setting
facts <- list(c(), 1:3)
units <- list(7, 4)
levels <- 3
etas <- list(1)
criteria <- c('I', 'D', 'A')
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


# test: MSOpt works ####
test_that("MSOpt works", {
  expect_equal(MSOpt(facts, units, levels, etas, criteria, model),
               list("facts" = list(c(), 1:3),
                    "nfacts" = 3,
                    "nstrat" = 2,
                    "units" = list(7, 4),
                    "runs" = 28,
                    "etas" = list(1),
                    "avlev" = list(c(-1, 0 , 1), c(-1, 0, 1), c(-1, 0, 1)),
                    "levs" = c(3, 3, 3),
                    "Vinv" = t(solve(diag(28) + 1 * kronecker(diag(7), matrix(1, 4, 4)))),
                    "model"  = 'quadratic',
                    "crit" = c('I', 'D', 'A'),
                    "ncrit" = 3,
                    "M" = Mq
               )
  )
})

# test: Score works ####

file_name <- here::here("tests/testthat/tests_data/ex_score_cs1.Rds")
# saveRDS(design_A, file = file_name)
example <- readRDS(file_name)
msopt <- MSOpt(facts, units, levels, etas, criteria, model)

test_that("Score works",{expect_equal(Score(msopt, example),
                                      c(0.3471977990, 0.1038650765, 0.1290310631),
                                      tolerance = 0.000000001
)
})



# test: MSSearch x I-optimality (CE) ####
set.seed(2)
msopt_I <- MSOpt(facts, units, levels, etas, c("I"), model)
#mssearch_I <- MSSearch(msopt_I, 1)
file_name <- here::here("tests/testthat/tests_data/mss_I_cs1.Rds")
#saveRDS(mssearch_I, file = file_name)
mss_I_cs1 <- readRDS(file_name)

set.seed(2)
test_that("MSSearch works", {
  expect_equal(MSSearch(msopt_I, 1),
               list("optsol" = mss_I_cs1$optsol,
                    "optsc" = mss_I_cs1$optsc,
                    "feval" = mss_I_cs1$feval,
                    "trend" = mss_I_cs1$trend
               )
  )
})

# test: MSSearch x I-optimality (CE) + parametri opz. ####
set.seed(2)
#mssearch_I2 <- MSSearch(msopt_I, 1, "Restarts", 50, "Start", example)
file_name <- here::here("tests/testthat/tests_data/mss_I2_cs1.Rds")
#saveRDS(mssearch_I2, file = file_name)
mss_I2_cs1 <- readRDS(file_name)

set.seed(2)
test_that("MSSearch works", {
  expect_equal(MSSearch(msopt_I, 1, "Restarts", 50, "Start", example),
               list("optsol" = mss_I2_cs1$optsol,
                    "optsc" = mss_I2_cs1$optsc,
                    "feval" = mss_I2_cs1$feval,
                    "trend" = mss_I2_cs1$trend
               )
  )
})

# test: MSSearch x I,D,A-optimality + parametri opz. ####
set.seed(2)
# mssearch3 <- MSSearch(msopt, c(2/4, 1/4, 1/4),
#                         "Restarts", 50,
#                         "Start", example,
#                         "Normalize", c(rep(0.5, 3), rep(1, 3)))
file_name <- here::here("tests/testthat/tests_data/mss_3_cs1.Rds")
#saveRDS(mssearch3, file = file_name)
mss_3_cs1 <- readRDS(file_name)

set.seed(2)
test_that("MSSearch works", {
  expect_equal(MSSearch(msopt, c(2/4, 1/4, 1/4),
                        "Restarts", 50,
                        "Start", example,
                        "Normalize", c(rep(0.5, 3), rep(1, 3))),
               list("optsol" = mss_3_cs1$optsol,
                    "optsc" = mss_3_cs1$optsc,
                    "feval" = mss_3_cs1$feval,
                    "trend" = mss_3_cs1$trend
               )
  )
})



# test: runTPLS ####

lCrit <- length(criteria)
iters <- 3 * lCrit
restarts <- 30
restInit <- 2

#tpls <- runTPLS(facts, units, criteria, model, iters,
#                "Restarts", restarts,
#                "RestInit", restInit,
#                "RngSeed", 4)
file_name <- here::here("tests/testthat/tests_data/tpls_cs1.Rds")
#saveRDS(tpls, file = file_name)
tpls_cs1 <- readRDS(file_name)

test_that("runTPLSearch works", {
  expect_equal(runTPLS(facts, units, criteria, model, iters,
                       "Etas", etas,
                       "Restarts", restarts,
                       "RestInit", restInit,
                       "RngSeed", 4),
               list("ar" = tpls_cs1$ar,
                    "stats" = tpls_cs1$stats,
                    "megaAR" = tpls_cs1$megaAR)
  )
}
)

