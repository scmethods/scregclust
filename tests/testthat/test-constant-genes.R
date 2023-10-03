test_that("constant genes are discarded correctly", {
  expression <- rbind(
    rep.int(1, 100),
    matrix(rnorm(500), nrow = 5),
    rep.int(0.5, 100),
    rnorm(100)
  )

  genesymbols <- c("T1", "T2", "T3", "T4", "T5", "T6", "R1", "R2")
  is_regulator <- c(0, 0, 0, 0, 0, 0, 1, 1)

  fit <- scregclust(
    expression, genesymbols, is_regulator, 0.1, 2, verbose = FALSE
  )

  expect_equal(
    fit$results[[1]]$genesymbols,
    c("T2", "T3", "T4", "T5", "T6", "R2")
  )
  expect_equal(
    fit$results[[1]]$is_regulator,
    c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  )
})
