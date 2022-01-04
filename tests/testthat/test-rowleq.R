a <- c(4, 5, 7, 6, 8)
b <- c(4, 5, 6, 7, 8)

test_that("rowleq works", {
  expect_equal(rowleq(a, b), FALSE)
})



