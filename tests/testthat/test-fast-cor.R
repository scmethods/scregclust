test_that("fast correlation computation", {
  pt <- 50
  pr <- 10
  n <- 200

  zt <- matrix(rnorm(pt * n), ncol = pt)
  zr <- matrix(rnorm(pr * n), ncol = pr)

  c_ref <- cor(zt, zr)
  c2 <- fast_cor(zt, zr)

  expect_equal(
    c_ref,
    c2
  )
})
