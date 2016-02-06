context("AAStringSet")

test_that(".readAAStringSet", {
  f <- file.path(system.file("extdata", package = "Pbase"),
                 "01_test_database.fasta")
  aa <- AAStringSet(c("P1" = "AKAKBKCKDE",
                      "P2" = "FKGHD",
                      "P3" = "JKKLMKNKDOPQRST"))
  df <- DataFrame(
    DB = Rle(factor(rep("td", 3))),
    AccessionNumber = paste0("P", 1:3),
    EntryName = paste0("P", 1:3, "_TEST"),
    IsoformName = Rle(rep(NA_character_, 3)),
    ProteinName = c("protein 1, length 10",
                    "protein 2, length 5",
                    "protein 3, length 15"),
    OrganismName = Rle(factor("machina arithmetica")),
    GeneName = Rle(factor(paste0("g", 1:3))),
    ProteinExistence = Rle(factor(rep(1, 3),
        labels = c("Evidence at protein level",
            "Evidence at transcript level",
            "Inferred from homology",
            "Predicted",
            "Uncertain"),
        levels = 1L:5L)),

    SequenceVersion = Rle(rep("1", 3)),
    Comment = Rle(c(rep(NA_character_, 3))),
    Filename = Rle(factor(rep(f, 3))))

  mcols(aa) <- df

  expect_identical(Pbase:::.readAAStringSet(f), aa)
  ## works with multiple files (and reordering them correctly)
  bb <- c(aa, aa)[c(1, 4, 2, 5, 3, 6)]
  pp <- Pbase:::.readAAStringSet(c(f, f))
  expect_identical(names(pp), names(bb))
  expect_identical(mcols(pp)$ProteinName, mcols(bb)$ProteinName)

  expect_error(Pbase:::.readAAStringSet("foobar"),
               "The file\\(s\\) .*foobar.* do\\(es\\) not exist!")
})
