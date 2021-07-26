set.seed(123)
options(digits = 10)

facts <- list(1, 2, 3)
units <- list(3, 3, 3)
levels <- 3
etas <- list(1, 1)
criteria <- "D"
model <- "interaction"

msopt <- MSOpt(facts, units, levels, etas, criteria, model)

sol <- matrix(c( -1, -1, -1,
                 -1, -1,  1,
                 -1, -1,  1,
                 -1,  1, -1,
                 -1,  1,  1,
                 -1,  1,  1,
                 -1,  1,  1,
                 -1,  1, -1,
                 -1,  1, -1,
                  1, -1, -1,
                  1, -1,  1,
                  1, -1, -1,
                  1,  1,  1,
                  1,  1, -1,
                  1,  1, -1,
                  1,  1, -1,
                  1,  1,  1,
                  1,  1,  1,
                 -1, -1, -1,
                 -1, -1,  1,
                 -1, -1,  1,
                 -1,  1, -1,
                 -1,  1,  1,
                 -1,  1,  1,
                 -1, -1, -1,
                 -1, -1,  1,
                 -1, -1, -1), ncol = 3, byrow = T)

# test_that("MSSearch works", {
#   expect_equal(MSSearch(msopt, 1, "Restarts", 100),
#                list("optsol" = sol,
#                     "optsc" = 0.1297854392,
#                     "feval" = 20420,
#                     "trend" = rep(0.1297854392, 100))
#                )
#   }
# )




