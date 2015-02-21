context("proteinCoverage-functions")

test_that(".proteinCoverage", {
    irl1 <- IRangesList(IRanges(start = 1, end = 5),
                        IRanges(start = 1, end = 16),
                        IRanges(start = 1, end = 10),
                        IRanges(start = 1, end = 20))
    irl2 <- IRangesList(IRanges(start = 1, end = 6),
                        IRanges(start = 13, end = 16),
                        IRanges(start = 15, end = 25))
    result <- c("P2" = 0, "P3" = 0.25, "P1" = 0.6, "P4" = NA)

    expect_error(Pbase:::.proteinCoverage(list(), irl1))
    expect_error(Pbase:::.proteinCoverage(irl2, list()))
    expect_error(Pbase:::.proteinCoverage(irl2, irl1),
                 "No names for .*pattern.* available!")
    names(irl2) <- rep("P3", 3)
    expect_error(Pbase:::.proteinCoverage(irl2, irl1),
                 "No duplicated names for .*pattern.* allowed!")
    names(irl2) <- c("P1", "P3", "P4")
    expect_error(Pbase:::.proteinCoverage(irl2, irl1),
                 "No names for .*subject.* available!")
    names(irl1) <- c("P3", "P3", "P1", "P4")
    expect_error(Pbase:::.proteinCoverage(irl2, irl1),
                 "No duplicated names for .*subject.* allowed!")
    names(irl1) <- c("P2", "P3", "P1", "P4")

    expect_equal(Pbase:::.proteinCoverage(irl2, irl1), result)
})

