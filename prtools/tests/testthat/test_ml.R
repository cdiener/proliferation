context("machine learning helpers")

test_that("metrics work", {
    x <- runif(100)
    y <- x*3 + rnorm(100)
    mod <- lm(y ~ x)
    y_hat <- predict(mod)
    m <- measures(y, y_hat)
    expect_equal(5, length(m))
    expect_true(all(m>0))
})

test_that("interactions work", {
    m <- matrix(1, nrow=20, ncol=20)
    colnames(m) <- 1:20
    ints <- inter(m)
    ecols <- 20*(20+1)/2
    expect_equal(c(20, ecols), dim(ints))
    expect_true(all.equal(unname(ints), matrix(1, 20, ecols)))
})
