test_that("Archive works", {
  expect_equal(Archive(3, 5),
               list("nsols" = 0,
                    "dim" = 3,
                    "scores" = matrix(0, 5, 3),
                    "solutions" = vector("list", 5)
               )
  )
})

ar <- Archive(3,5)

test_that("Resize works", {
  expect_equal(Resize(ar),
               list("nsols" = 0,
                    "dim" = 3,
                    "scores" = matrix(0, 10, 3),
                    "solutions" = vector("list", 10)
               )
  )
})

ar <- Resize(ar)
sol <- "sol 1"
score <- c(324, 55, 6)

test_that("Add works", {
  expect_equal(Add(ar, sol, score),
               list("nsols" = 1,
                    "dim" = 3,
                    "scores" = matrix(c(324, 55, 6, rep(0, 27)), 10, 3,
                                      byrow = T),
                    "solutions" = c(sol, vector("list", 9))
               )
  )
}
)

data <- matrix(c(3, 3, 2, 2), ncol = 2, byrow = 2)

test_that("FixRepmat works", {
  expect_equal(FixRepmat(data, 2, 3),
               matrix(c(3, 3, 3, 3, 3, 3,
                        2, 2, 2, 2, 2, 2,
                        3, 3, 3, 3, 3, 3,
                        2, 2, 2, 2, 2, 2),
                      ncol = 6, byrow = T
               )
  )
}
)

data <- c(1, 2, 3)
test_that("FixRepmat works", {
  expect_equal(FixRepmat(data, 5, 1),
               matrix(rep(c(1, 2, 3), 5),
                      ncol = 3, byrow = T
               )
  )
})

test_that("FixRepmat works", {
  expect_equal(FixRepmat(8, 1, 3), matrix(c(8,8,8), ncol = 3))
})

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

test_that("RemoveDominated works", {
  expect_equal(RemoveDominated(ar),
               list("nsols" = 6,
                    "dim" = 3,
                    "scores" = matrix(c( 1.0, 2, 3,
                                         0.5, 3, 4,
                                         0.5, 3, 3,
                                         2.0, 1, 1,
                                         2.0, 3, 3,
                                         2.0, 3, 3),
                                      ncol = 3, byrow = T
                    ),
                    "solutions" = list("sol 2",
                                       "sol 4", "sol 6",
                                       "sol 8", "sol 9",
                                       "sol 10")
               )
  )
})

ar = RemoveDominated(ar)

test_that("RemoveDuplicates works", {
  expect_equal(RemoveDuplicates(ar),
               list("nsols" = 5,
                    "dim" = 3,
                    "scores" = matrix(c( 1.0, 2, 3,
                                         0.5, 3, 4,
                                         0.5, 3, 3,
                                         2.0, 1, 1,
                                         2.0, 3, 3),
                                      ncol = 3, byrow = T),
                    "solutions" = list("sol 2", "sol 4",
                                       "sol 6", "sol 8",
                                       "sol 9"))
  )
})

# ar = RemoveDuplicates(ar)
# ar = RemoveDominated(ar)

ar1 = Archive(3,3)
ar1 = Resize(ar1)

ar1 = Add(ar1, "sol 1", c(4, 5, 6))
ar1 = Add(ar1, "sol 2", c(1, 2, 3))
ar1 = Add(ar1, "sol 3", c(1, 1, 6))
ar1 = Add(ar1, "sol 4", c(0.5, 3, 4))

test_that("Trim works", {
  expect_equal(Trim(ar1),
               list("nsols" = 4,
                    "dim" = 3,
                    "scores" = matrix(c(4, 5, 6,
                                        1, 2, 3,
                                        1, 1, 6,
                                        0.5, 3, 4),
                                      4, 3, byrow = T
                    ),
                    "solutions" = list("sol 1", "sol 2",
                                       "sol 3", "sol 4"
                    )
               )
  )
})
