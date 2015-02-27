context("Proteins-functions")

test_that("proteinCoverage", {
    # build fake Proteins
    p <- Proteins()
    p@aa <- AAStringSet(c(P1="ABCDEF", P2="ABCD", P3="HIJKL"))
    p@pranges <- IRangesList(P1=IRanges(1, 6), P2=IRanges(1, 3), P3=IRanges())
    mcols(p@aa) <- DataFrame(AccessionNumber=c("P1", "P2", "P3"))

    expect_equal(acols(proteinCoverage(p))$Coverage, c(P1=1, P2=0.75, P3=0))
})

