test_that("Set parameter space creates all possible combinations", {
    a = b <- round(seq(1, 4, 1), 3)
    alpha = gamma <- round(seq(1, 3, 1))
    expect_equal(nrow(set_parameter_space(a, b, alpha, gamma)),
    round(sum(seq(1, length(a)-1, 1)), 1)*length(alpha)*length(gamma))
})
