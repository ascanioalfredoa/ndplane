test_that("trapz calculates area correctly", {
  x <- 0:10
  y <- rep(1, 11)
  expect_equal(trapz(x, y), 10)

  x2 <- c(0, 1, 2)
  y2 <- c(0, 1, 0)
  expect_equal(trapz(x2, y2), 1)
})

test_that("beta_overlap merges and calculates min y", {
  spa <- data.frame(x = 0:10, y = rep(1, 11))
  spb <- data.frame(x = 0:10, y = seq(0, 1, 0.1))
  ov <- beta_overlap(spa, spb)
  expect_equal(ov$y, spb$y)
})

test_that("niche indices are within [0, 1]", {
  spa <- betaPDF(a = 0, b = 1, alpha = 2, gamma = 2)
  spb <- betaPDF(a = 0.5, b = 1.5, alpha = 2, gamma = 2)
  ov <- beta_overlap(spa, spb)

  diss <- niche_diss(spa, spb, ov)
  excl <- niche_excl(spa, spb)

  expect_true(diss >= 0 && diss <= 1)
  expect_true(excl >= 0 && excl <= 1)
})
