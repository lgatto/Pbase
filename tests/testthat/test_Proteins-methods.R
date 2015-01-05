context("Proteins-methods")

f <- file.path(system.file("extdata", package = "Pbase"),
               "01_test_database.fasta")

p <- Proteins(f)

test_that("cleave", {
    pc <- cleave(p)
    expect_identical(length(pc@pranges), 3L)
    expect_identical(pc@pranges@unlistData@elementMetadata[, "MissedCleavages"],
                     Rle(0, 12))
    pc <- cleave(p, missedCleavages = 2)
    expect_identical(pc@pranges@unlistData@elementMetadata[, "MissedCleavages"],
                     Rle(2, 7))
})

test_that("isCleaved", {
    expect_false(isCleaved(p))
    expect_false(isCleaved(cleave(p, missedCleavages = 1), missedCleavages = 2))
    expect_true(isCleaved(cleave(p)))
    expect_true(isCleaved(cleave(p, missedCleavages = 2), missedCleavages = 2))
})

