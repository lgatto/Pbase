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

test_that("addPeptideFragments", {
    fragments <- tempfile(fileext=".fasta")
    on.exit(unlink(fragments))

    irl <- IRangesList(P1=IRanges(c(8, 1), c(10, 4)),
                       P2=IRanges(),
                       P3=IRanges(c(2, 10), c(4, 15)))

    writeLines(c(">td|P1|PEP1 peptide 1, length 3 OS=machina arithmetica GN=g1 PE=1 SV=1",
                 "KDE",
                 ">td|P1|PEP2 peptide 2, length 4 OS=machina arithmetica GN=g1 PE=1 SV=1",
                 "AKAK",
                 ">td|P1|PEP3 peptide 3, length 4 OS=machina arithmetica GN=g1 PE=1 SV=1",
                 "ZZZ",
                 ">td|P3|PEP4 peptide 4, length 3 OS=machina arithmetica GN=g3 PE=1 SV=1",
                 "KKL",
                 ">td|P4|PEP5 peptide 5, length 3 OS=machina arithmetica GN=g4 PE=1 SV=1",
                 "ABC",
                 ">td|P3|PEP6 peptide 6, length 6 OS=machina arithmetica GN=g3 PE=1 SV=1",
                 "OPQRST"), fragments)

    df <- DataFrame(
        DB = Rle(factor(rep("td", 4))),
        AccessionNumber = paste0("P", rep(c(1, 3), each=2)),
        EntryName = paste0("PEP", c(1:2, 4, 6)),
        IsoformName = Rle(rep(NA_character_, 4)),
        ProteinName = c("peptide 1, length 3",
                        "peptide 2, length 4",
                        "peptide 4, length 3",
                        "peptide 6, length 6"),
        OrganismName = Rle(factor("machina arithmetica")),
        GeneName = Rle(factor(rep(c(1, 2), each=2),
                              labels = paste0("g", c(1, 3:4)),
                              levels = 1:3)),
        ProteinExistence = Rle(factor(rep(1, 4),
            labels = c("Evidence at protein level",
                "Evidence at transcript level",
                "Inferred from homology",
                "Predicted",
                "Uncertain"),
            levels = 1L:5L)),
        SequenceVersion = Rle(rep("1", 4)),
        Comment = Rle(c(rep(NA_character_, 4))),
        Filename = Rle(factor(rep(fragments, 4))),
        PeptideIndex = Rle(c(1, 2, 4, 5)),
        ProteinIndex = Rle(rep(c(1, 3), each=2)))

    pr <- pranges(addPeptideFragments(p, fragments))
    expect_true(length(pr) == 2)
    expect_true(all(unlist(pr) == unlist(irl)))
    expect_equal(mcols(pr[[1]]), df[1:2, ])
    expect_equal(mcols(pr[[2]]), df[3:4, ])

    pr <- pranges(addPeptideFragments(p, fragments, rmEmptyRanges=FALSE))
    expect_true(length(pr) == 3)
    expect_true(all(unlist(pr) == unlist(irl)))
    expect_equal(mcols(pr[[1]]), df[1:2, ])
    expect_equal(mcols(pr[[2]]), df[0, ])
    expect_equal(mcols(pr[[3]]), df[3:4, ])

    expect_error(addPeptideFragments(p, "foobar"),
                 "The file\\(s\\) .*foobar.* do\\(es\\) not exist!")
})

test_that("pranges replacement", {
    expect_error(pranges(p) <- 1:3,
                 "unable to find an inherited method for function .*pranges<-.* for signature .*Proteins.*, .*integer.*")
    expect_error(pranges(p) <- IRangesList(),
                 "Length of replacement pranges differs from current ones.")
    expect_error(pranges(p) <- IRangesList(A=IRanges(1, 2), B=IRanges(1, 2), C=IRanges(1, 2)),
                 "Names of replacement pranges differ from current ones.")

    pm <- p
    irl <- IRangesList(P1=IRanges(1, 2),
                       P2=IRanges(2, 3),
                       P3=IRanges(3, 4))
    pranges(pm) <- irl
    expect_equal(pranges(pm), irl)
    expect_error(pranges(p) <- irl[3:1],
                 "Names of replacement pranges differ from current ones.")

    pc <- cleave(p)
    pranges(pm) <- pranges(pc)
    expect_equal(pranges(pm), pranges(pc))
    l <- LogicalList(c(TRUE, FALSE, FALSE, TRUE),
                     c(TRUE, FALSE),
                     c(rep(TRUE, 3), rep(FALSE, 3)))
    pranges(pm) <- pranges(pm)[l]
    expect_equal(pranges(pm), pranges(pc)[NumericList(c(1, 4), c(1), 1:3)])
})

test_that("acols replacement", {
    expect_error(acols(p) <- 1:3,
                 "unable to find an inherited method for function .*acols<-.* for signature .*Proteins.*, .*integer.*")
    expect_error(acols(p) <- DataFrame(),
                 "Number of rows of replacement acols differ from current ones.")

    pm <- p
    ac <- DataFrame(A=1:3, B=1:3, row.names = c("P1", "P2", "P3"))
    acols(pm) <- ac
    rownames(pm@aa@elementMetadata) <- c("P1", "P2", "P3")
    expect_equal(acols(pm), ac)
    expect_error(acols(pm) <- ac[3:1,],
                 "Row names of replacement acols differ from current ones.")
})

## Unit test for issue #27; thanks to Johannes Rainer (@jotsetung) for
## reporting and fixing
test_that("pmetadata", {
    ## Create the pranges: have a IRange for the 1st and 3rd.
    ir <- IRanges(start = c(3, 5), end = c(10, 15))
    mcols(ir) <- DataFrame(AccessionNumber = c("P1", "P3"), OtherMcol = c(1, 3))
    irL <- split(ir, f = mcols(ir)$AccessionNumber)
    ## Add the empty one:
    emptyIr <- IRanges()
    mcols(emptyIr) <- DataFrame(
        matrix(ncol = 2, nrow = 0,
               dimnames = list(rownames = character(),
                               colnames = c("AccessionNumber", "OtherMcol"))))
    ## Create the IRangesList
    irL <- c(irL[1], IRangesList(emptyIr), irL[2])
    names(irL) <- c("P1", "P2", "P3")
    ## Add the IRangesList to the Proteins object
    pranges(p) <- irL

    ## We have 3 sequences, thus we should expect 3 elements:
    expect_true(length(pcols(p)) == length(p))
    expect_true(length(p@pranges) == length(p))
    ## Names and order should match
    expect_identical(names(pcols(p)), seqnames(p))
    expect_identical(names(p@pranges), seqnames(p))
    ## The length of the pcols and pranges have to match
    expect_identical(lengths(pranges(p)), lengths(pcols(p)))
    ## nrow of the pcols should be 1, 0, 1:
    expect_equal(lengths(pcols(p)), setNames(c(1, 0, 1), c("P1", "P2", "P3")))
    expect_equal(elementNROWS(pranges(p)), setNames(c(1, 0, 1), c("P1", "P2", "P3")))
})
