context("utils")

test_that(".flatIRangesList", {
  expect_error(Pbase:::.flatIRangesList(list()))

  irl <- IRangesList(IRanges(start = c(1, 5), end = c(3, 10)),
                     IRanges(start = c(7, 10), end = c(9, 15)))
  result <- IRanges(start = c(1, 5, 7, 10), end = c(3, 10, 9, 15))
  resultShift <- IRanges(start = c(1, 5, 17, 20), end = c(3, 10, 19, 25))

  expect_equal(Pbase:::.flatIRangesList(irl), result)
  expect_equal(Pbase:::.flatIRangesList(irl, shift = TRUE), resultShift)

})

test_that(".splitIRanges", {
  expect_error(Pbase:::.splitIRanges(list()))

  ir <- IRanges(start = c(1, 5, 17, 20), end = c(3, 10, 19, 25))

  result <- IRangesList(IRanges(start = 1, end = 3),
                        IRanges(start = 5, end = 10),
                        IRanges(start = 17, end = 19),
                        IRanges(start = 20, end = 25))

  resultSplit <- IRangesList(IRanges(start = c(1, 5), end = c(3, 10)),
                             IRanges(start = c(17, 20), end = c(19, 25)))
  names(resultSplit) <- as.character(1:2)

  resultUnshift <- IRangesList(IRanges(start = c(1, 5), end = c(3, 10)),
                               IRanges(start = c(7, 10), end = c(9, 15)))

  expect_equal(Pbase:::.splitIRanges(ir, USE.NAMES = FALSE), result)
  expect_equal(Pbase:::.splitIRanges(ir), setNames(result, as.character(1:4)))
  expect_equal(Pbase:::.splitIRanges(ir, f = c(1, 1, 2, 2)), resultSplit)
  expect_equal(Pbase:::.splitIRanges(ir, f = c(1, 1, 2, 2), unshift = TRUE,
                                     USE.NAMES = FALSE), resultUnshift)
})

