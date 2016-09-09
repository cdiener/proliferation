context("analysis tools")

test_that("tissue lfcs can be calculated", {
    v <- matrix(1, nrow=10, ncol=10)
    colnames(v) <- paste0("v", 1:10)
    rownames(v) <- paste0("s", 1:10)
    map <- rep(1:2, each=5)
    names(map) <- rownames(v)
    lf <- panel_lfc(v, map)
    expect_equal(c(20, 3), dim(lf))
    extra <- data.frame(info=paste0("info", 1:10))
    lf <- panel_lfc(v, map, extra)
    expect_equal(c(20, 4), dim(lf))
})

test_that("ES works as expected", {
    pws <- c("a", rep("b", 99))
    w <- rep(1, 100)
    expect_true(ES("a", w, pws) > 0)
    expect_true(ES("a", w, rev(pws)) < 0)
    expect_equal(2, length(ES("a", w, pws, both=TRUE)))
})

test_that("NES works as expected", {
    pws <- c("a", rep("b", 99))
    w <- rep(1, 100)
    nes <- NES("a", w, pws)
    expect_true(nes[1] > 0)
    expect_true(nes[2] < 0.05)
    nes <- NES("a", w, rev(pws))
    expect_true(nes[1] < 0)
    expect_true(nes[2] < 0.05)
})

test_that("shorten works", {
    short <- "abc"
    long <- paste(letters, collapse=" ")
    expect_equal(3, nchar(shorten(short, 5)))
    expect_equal(12, nchar(shorten(long, 10)))
})
