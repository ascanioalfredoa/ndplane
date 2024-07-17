test_that("betafuncion_max() produces an R expression", {
  expect_true(is.expression(betafunction_max()))
})

test_that("betaPDF produces a data frame", {
    expect_true(is.data.frame(betaPDF()))
})

