### test MSSearch multicrit - model = "quadratic" ####

# instance 1
set.seed(123)
options(digits = 10)

facts <- list(1, 2:5)
units <- list(6, 5)
levels <- 3
etas <- list(1)
criteria <- c('I', 'Id', 'D', 'A', 'Ds', 'As')
model <- "quadratic"

msopt <- MSOpt(facts, units, levels, etas, criteria, model)
load("C:\\Users\\Francesca\\Desktop\\Rtesi\\multiDoE\\mssearchM_i1.RData")

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt, rep(1/6, 6), "Restarts", 100),
               list("optsol" = mssearchM$optsol,
                    "optsc" = mssearchM$optsc,
                    "feval" = mssearchM$feval,
                    "trend" = mssearchM$trend
               )
  )
}
)

### test MSSearch multicrit - model = "interaction" ####
set.seed(12)
msopt <- MSOpt(facts, units, levels, etas, criteria, "interaction")
load("C:\\Users\\Francesca\\Desktop\\Rtesi\\multiDoE\\mss_int.RData")

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt, rep(1/6, 6), "Restarts", 100),
               list("optsol" = mss_int$optsol,
                    "optsc" = mss_int$optsc,
                    "feval" = mss_int$feval,
                    "trend" = mss_int$trend
               )
  )
}
)

### test MSSearch multicrit - model = "interaction", initial sol ####
set.seed(34)
# example
ex <- matrix(c( 0,  1,  0,  0,  0,
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
# solution
solution <- matrix(c(-1, 1, 1,-1,-1,
                      -1, 1,-1, 1, 1,
                      -1, 1,-1,-1,-1,
                      -1,-1,-1, 1,-1,
                      -1,-1, 1, 1, 1,
                      1, 1,-1,-1, 1,
                      1, 1, 1, 1, 1,
                      1,-1, 1,-1, 1,
                      1,-1, 1, 1,-1,
                      1, 1, 0,-1,-1,
                      -1,-1, 1,-1, 0,
                      -1, 1,-1,-1, 1,
                      -1,-1,-1, 1, 1,
                      -1, 1, 1, 1, 1,
                      -1, 1, 1, 1,-1,
                      1,-1,-1, 1,-1,
                      1,-1,-1,-1,-1,
                      1, 1, 1, 1,-1,
                      1,-1, 1, 1, 1,
                      1, 1,-1, 1, 1,
                      -1,-1,-1,-1, 1,
                      -1, 1,-1, 1,-1,
                      -1, 1, 1,-1, 1,
                      -1,-1,-1,-1,-1,
                      -1,-1, 1, 1,-1,
                      1, 1, 1,-1, 1,
                      1,-1, 1,-1,-1,
                      1, 1,-1, 1,-1,
                      1,-1,-1,-1, 1,
                      1,-1,-1, 1, 1),
                      ncol = 5, byrow = T)

test_that("MSSearch works", {
  expect_equal(MSSearch(msopt, rep(1/6, 6), "Restarts", 100,
                        "Start", ex),
               list("optsol" = solution,
                    "optsc" = c(0.359627855309882, 0.158847044392840,
                                0.045729043820545, 0.058181333743426,
                                0.041444756799279, 0.04867470193185),
                    "feval" = 25834,
                    "trend" = rep(0.118750789332971, 100)
               ),
               tolerance = 0.000000000001
  )
}
)

