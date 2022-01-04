# test: X matrix ####
X = matrix(c(3, 3, 3, 3, 2, 3, 2, 2), nrow = 2, byrow = T)
test_that("colprod works", {
  expect_equal(colprod(X),
               matrix(c(9, 9, 9, 9, 9, 9, 6, 4, 4, 6, 6, 4),
                      nrow = 2, byrow = T)
  )
})

# test: X vector ####
X = c(2, 2, 2, 5, 7, 9)
test_that("colprod works", {
  expect_equal(colprod(X),
               c(4, 4, 10, 14, 18, 4, 10, 14, 18, 10, 14, 18, 35, 45, 63))
})

# test: X scalar ####
X = 3
test_that("colprod works", {
  expect_equal(colprod(X), NULL)
})

