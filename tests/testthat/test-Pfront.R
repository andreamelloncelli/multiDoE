# TEST PFront ####
ar = Archive(3,5)
ar = Resize(ar)
ar = Add(ar, "sol 1", c(4, 5, 6))
ar = Add(ar, "sol 2", c(1, 2, 3))
ar = Add(ar, "sol 3", c(3, 3, 6))
ar = Add(ar, "sol 4", c(0.5, 3, 4))
ar = Add(ar, "sol 5", c(3, 2, 6))
ar = Add(ar, "sol 6", c(0.5, 3, 3))
ar = Add(ar, "sol 7", c(3, 3, 6))
ar = Add(ar, "sol 8", c(2, 1, 1))
ar = Add(ar, "sol 9", c(2, 3, 3))
ar = Add(ar, "sol 10", c(2, 3, 3))

#ar = RemoveDominated(ar)

pf = PFront(ar)

addpf = Add_PF(pf, ar$nsols)

test_that("PFront works",
          {expect_equal(PFront(ar), list("arch" = ar,
                                         "gaps" = list(),
                                         "scmax" = c(2, 1, 1),
                                         "scmin" = c(2, 1, 1),
                                         "ptrs" = 8
                            )
                        )
          }
)




# ####


SetMinMaxSc <- function(pf, scmax, scmin) {
  pf$scmax <- scmax
  pf$scmin <- scmin
  return(pf)
} # OK



# TEST dominatedHypervolume ####
options(digits = 10)
n = 200
theta = seq( - 3 * pi / 2, - pi, length.out = n)
FP = cbind(1 + cos(theta), 1 - sin(theta))

test_that("dominatedHypervolume works",
          {expect_equal(dominatedHypervolume(FP, c(1, 1)) , 0.783416635050702)
          }
)


