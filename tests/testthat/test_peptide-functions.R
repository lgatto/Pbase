context("peptides-functions")

test_that(".calculateMolecularWeight", {
    x <- c(ACE = "ACE", ACDEFGHIKLMN = "ACDEFGHIKLMN")
    result <- c(ACE = 321.1, ACDEFGHIKLMN = 1376.6)
    expect_equal(round(Pbase:::.calculateMolecularWeight(x), 1), result)
})

test_that(".isValidPeptide", {
    x <- c(ACE = "ACE", ACDEFGHIKLMN = "ACDEFGHIKLMN")
    expect_equal(Pbase:::.isValidPeptide(x), c(ACE = TRUE, ACDEFGHIKLMN = TRUE))
    expect_equal(unname(Pbase:::.isValidPeptide(x, len = c(4, 6))),
                       c(FALSE, FALSE))
    expect_equal(unname(Pbase:::.isValidPeptide(x, len = c(4, 20))),
                 c(FALSE, TRUE))
    expect_equal(unname(Pbase:::.isValidPeptide(x, mass = c(0, 10),
                                                len = c(4, 20))),
                 c(FALSE, FALSE))
    expect_equal(unname(Pbase:::.isValidPeptide(x, mass = c(300, 1400),
                                                len = c(1, 20))), c(TRUE, TRUE))
})

test_that(".peptidePosition", {
    pattern <- c(P1 = "ACE", P2 = "IL", P3 = "KR")
    subject <- c(P1 = "ACEDE", P3 = "JKLMN", P2 = "FGHIL")
    result <- IRangesList(c(P1 = IRanges(1, 3),
                            P2 = IRanges(4, 5)))
    expect_error(Pbase:::.peptidePosition(pattern,
                                          setNames(subject, c(rep("P1", 3)))),
                 "No duplicated names for .*subject.* allowed!")
    expect_equal(Pbase:::.peptidePosition(pattern, subject), result)
})

