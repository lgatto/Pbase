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
    ## one peptide per protein
    pattern <- c(P1 = "ACE", P2 = "IL", P3 = "KR")
    subject <- c(P1 = "ACEDE", P3 = "JKLMN", P2 = "FGHIL")
    result <- IRangesList(c(P1 = IRanges(1, 3),
                            P2 = IRanges(4, 5)))
    result <- IRanges(start = c(1,   # P1
                                4),  # P2
                      end = c(3,     # P1
                              5))    # P2
    mcols(result) <- DataFrame(PeptideIndex = Rle(c(1, 2)),
                               ProteinIndex = Rle(c(1, 3)))
    result <- split(result, f = c("P1", "P2"))

    expect_error(Pbase:::.peptidePosition(unname(pattern), subject),
                 "No names for .*pattern.* available!")
    expect_error(Pbase:::.peptidePosition(pattern, unname(subject)),
                 "No names for .*subject.* available!")
    expect_error(Pbase:::.peptidePosition(pattern,
                                          setNames(subject, c(rep("P1", 3)))),
                 "No duplicated names for .*subject.* allowed!")

    expect_equal(Pbase:::.peptidePosition(pattern, subject), result)

    ## multiple peptides per protein
    ## see https://github.com/ComputationalProteomicsUnit/Pbase/issues/5
    pattern <- c(P1 = "DE", P2 = "IL", P3 = "KR", P1 = "ABC", P3 = "LMN")
    pattern2 <- AAStringSetList(P1 = AAStringSet(c("DE", "ABC")),
                                P2 = AAStringSet("IL"),
                                P3 = AAStringSet(c("KR", "LMN")))
    subject <- c(P1 = "ABCDE", P3 = "JKLMN", P2 = "FGHIL")
    result <- IRanges(start = c(4, 1,   # P1
                                4,      # P2
                                3),     # P3
                      end = c(5, 3,     # P1
                              5,        # P2
                              5))       # P3
    result2 <- result
    mcols(result) <- DataFrame(PeptideIndex = Rle(c(1, 4, 2, 5)),
                               ProteinIndex = Rle(c(1, 1, 3, 2)))
    mcols(result2) <- DataFrame(PeptideIndex = Rle(c(1, 2, 3, 5)),
                                ProteinIndex = Rle(c(1, 1, 3, 2)))
    result <- split(result, f = c("P1", "P1", "P2", "P3"))
    result2 <- split(result2, f = c("P1", "P1", "P2", "P3"))
    expect_equal(Pbase:::.peptidePosition(pattern, subject), result)
    expect_equal(Pbase:::.peptidePosition(pattern2, subject), result2)
})

