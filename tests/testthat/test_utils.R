context("utils")

test_that(".addColumn", {
    d <- data.frame(a = 1:3, b = 4:6)
    result <- data.frame(a = 1:3, b = 4:6, c = 7:9)
    resultForced <- data.frame(a = 1:3, b = 7:9)

    expect_error(Pbase:::.addColumn(list()))
    expect_error(Pbase:::.addColumn(d))
    expect_error(Pbase:::.addColumn(d, character()))
    expect_error(Pbase:::.addColumn(d, ""))
    expect_error(Pbase:::.addColumn(d, "b", integer()))
    expect_error(Pbase:::.addColumn(d, "b", 7:9),
                 "The column .*b.* already exists.")

    expect_equal(Pbase:::.addColumn(d, "c", 7:9), result)
    expect_equal(Pbase:::.addColumn(d, "b", 7:9, force = TRUE), resultForced)
})

test_that(".isInRange", {
    x <- 1:10
    result <- rep(c(FALSE, TRUE, FALSE), c(2, 4, 4))
    expect_equal(Pbase:::.isInRange(x, 3, 6), result)
    expect_equal(Pbase:::.isInRange(x, 1, 20), rep(TRUE, 10))
    expect_equal(Pbase:::.isInRange(x, 11, 20), rep(FALSE, 10))
})

test_that(".flatIRangesList", {
    expect_error(Pbase:::.flatIRangesList(list()))

    irl <- IRangesList(IRanges(start = c(1, 5), end = c(3, 10)),
                       IRanges(start = c(7, 10), end = c(9, 15)))
    result <- IRanges(start = c(1, 5, 7, 10), end = c(3, 10, 9, 15))
    resultShift <- IRanges(start = c(1, 5, 17, 20), end = c(3, 10, 19, 25))
    resultShiftBy <- IRanges(start = c(2, 6, 27, 30), end = c(4, 11, 29, 35))

    expect_equal(Pbase:::.flatIRangesList(irl), result)
    expect_equal(Pbase:::.flatIRangesList(irl, shift = TRUE), resultShift)
    expect_equal(Pbase:::.flatIRangesList(irl, shift = TRUE,
                                          shiftBy = c(1, 20)), resultShiftBy)

})

test_that(".splitIRanges", {
    expect_error(Pbase:::.splitIRanges(list()))

    ir <- IRanges(start = c(1, 5, 17, 20), end = c(3, 10, 19, 25))

    result <- IRangesList(IRanges(start = 1, end = 3),
                          IRanges(start = 5, end = 10),
                          IRanges(start = 17, end = 19),
                          IRanges(start = 20, end = 25))
    names(result) <- as.character(1:4)

    resultSplit <- IRangesList(IRanges(start = c(1, 5), end = c(3, 10)),
                               IRanges(start = c(17, 20), end = c(19, 25)))
    names(resultSplit) <- as.character(1:2)

    resultUnshift <- IRangesList(IRanges(start = c(1, 5), end = c(3, 10)),
                                 IRanges(start = c(7, 10), end = c(9, 15)))
    names(resultUnshift) <- as.character(1:2)

    expect_equal(Pbase:::.splitIRanges(ir), result)
    expect_equal(Pbase:::.splitIRanges(ir, f = c(1, 1, 2, 2)), resultSplit)
    expect_equal(Pbase:::.splitIRanges(ir, f = c(1, 1, 2, 2), unshift = TRUE),
                                       resultUnshift)
})

test_that(".setNames2", {
    x <- 1:3
    nm1 <- LETTERS[1:3]
    nm2 <- setNames(LETTERS[1:3], LETTERS[4:6])
    expect_equal(Pbase:::.setNames2(x, nm1), setNames(x, nm1))
    expect_equal(Pbase:::.setNames2(x, nm2), setNames(x, names(nm2)))
})

test_that(".singleAA", {
    x1 <- "ABC"
    x2 <- c(a = "ABC", b = "DEF")
    x3 <- AAString(x1)
    x4 <- AAStringSet(x2)
    result1 <- list(LETTERS[1:3])
    result2 <- list(a = LETTERS[1:3], b = LETTERS[4:6])
    expect_equal(Pbase:::.singleAA(x1), result1)
    expect_equal(Pbase:::.singleAA(x2), result2)
    expect_equal(Pbase:::.singleAA(x3), result1)
    expect_equal(Pbase:::.singleAA(x4), result2)
})

test_that(".singular", {
    x1 <- c(1, 2, 2, 3, 1, 5, 4, 3, 2, 1)
    x2 <- rep(1, 10)
    x3 <- 1:5
    expect_equal(Pbase:::.singular(x1), c(5, 4))
    expect_equal(Pbase:::.singular(x2), double())
    expect_equal(Pbase:::.singular(x3), 1:5)
})

