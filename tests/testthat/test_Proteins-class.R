context("Proteins-class")

test_that("validity", {
  expect_error(new("Proteins", seq=AAStringSet("ABC")),
               paste0("Number of rows in the metadata and ",
                      "the length of the sequences do not match!"))
})

test_that("fasta constructor", {
  expect_error(Proteins("foobar.txt"),
               "The file.s. .*foobar.txt.* do.es. not exist!")
               #paste0("The file(s) ", sQuote("foobar.txt"),
               #       " do(es) not exist!"))

  f <- file.path(system.file("extdata", package = "Pbase"),
                 "01_test_database.fasta")
  aa <- AAStringSet(c("sp|P1|protein 1, length 10"="AKAKBKCKDE",
                      "sp|P2|protein 2, length 5"="FKGHD",
                      "sp|P3|protein 3, length 15"="JKKLMKNKDOPQRST"))
  p <- Proteins(f)
  p2 <- p
  p2@seq <- p@seq[1:2]
  p2@mcols <- p@mcols[1:2,, drop=FALSE]

  expect_identical(length(p), 3L)
  expect_identical(nrow(mcols(p)), 3L)
  expect_identical(as.integer(mcols(p)[1L, ]), 1L)
  expect_identical(basename(levels(mcols(p)[1L, ])), "01_test_database.fasta")

  ## could not be compare directly because testthat throws an error:
  ##  "cannot unclass an external pointer"
  expect_identical(as.character(p@seq), as.character(aa))
  expect_identical(as.character(p[[1]]), as.character(aa[[1]]))
  expect_identical(as.character(p[1:2]@seq), as.character(p2@seq))
  expect_identical(p[1:2]@mcols, p2@mcols)

  expect_true(is.character(metadata(p)$created))
  expect_identical(nchar(metadata(p)$created,), 24L)
  expect_true(grepl(paste0("[A-z]{3} [A-z]{3} +[0-9]{1,2} ",
                           "[0-9]{2}:[0-9]{2}:[0-9]{2} [0-9]{4}"),
                    metadata(p)$created))
})

