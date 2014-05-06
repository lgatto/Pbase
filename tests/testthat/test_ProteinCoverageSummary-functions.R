context("ProteinCoverageSummary-functions")

test_that(".calculateProteinCoverage", {
  irl1 <- IRangesList(IRanges(start = 1, end = 5),
                      IRanges(start = 1, end = 15),
                      IRanges(start = 1, end = 10))
  irl2 <- IRangesList(IRanges(start = 1, end = 6),
                      IRanges(start = 11, end = 15))
  result <- LogicalList(P2 = logical(5),
                        P3 = rep(c(FALSE, TRUE), c(10, 5)),
                        P1 = rep(c(TRUE, FALSE), c(6, 4)))
  expect_error(Pbase:::.calculateProteinCoverage(list(), irl1))
  expect_error(Pbase:::.calculateProteinCoverage(irl2, list()))
  expect_error(Pbase:::.calculateProteinCoverage(irl2, irl1),
               "No names for .*pattern.* available!")
  names(irl2) <- c("P3", "P3")
  expect_error(Pbase:::.calculateProteinCoverage(irl2, irl1),
               "No duplicated names for .*pattern.* allowed!")
  names(irl2) <- c("P1", "P3")
  expect_error(Pbase:::.calculateProteinCoverage(irl2, irl1),
               "No names for .*subject.* available!")
  names(irl1) <- c("P3", "P3", "P1")
  expect_error(Pbase:::.calculateProteinCoverage(irl2, irl1),
               "No duplicated names for .*subject.* allowed!")
  names(irl1) <- c("P2", "P3", "P1")

  expect_identical(Pbase:::.calculateProteinCoverage(irl2, irl1),
                   result)
})

