context("aa-functions")

test_that(".calculateMolecularWeight", {
    x <- c(ACE = "ACE", ACDEFGHIKLMN="ACDEFGHIKLMN")
    result <- c(ACE = 321.1, ACDEFGHIKLMN=1376.6)
    expect_equal(round(Pbase:::.calculateMolecularWeight(x), 1), result)
})

