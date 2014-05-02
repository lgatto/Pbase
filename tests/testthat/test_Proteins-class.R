context("Proteins-class")

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
  p2@aa <- p@aa[1:2]
  p2@aa@elementMetadata <- p@aa@elementMetadata[1:2,, drop=FALSE]

  expect_identical(length(p), 3L)
  expect_identical(nrow(ametadata(p)), 3L)
  expect_identical(as.integer(ametadata(p)[1L, ]), 1L)
  expect_identical(basename(levels(ametadata(p)[1L, ])),
                   "01_test_database.fasta")

  ## could not be compare directly because testthat throws an error:
  ##  "cannot unclass an external pointer"
  expect_identical(as.character(p@aa), as.character(aa))
  expect_identical(as.character(p[[1]]), as.character(aa[[1]]))
  expect_identical(as.character(p[1:2]@aa), as.character(p2@aa))

  expect_identical(ametadata(p[1:2]), ametadata(p2))

  expect_true(is.character(metadata(p)$created))
  expect_identical(nchar(metadata(p)$created,), 24L)
  expect_true(grepl(paste0("[A-z]{3} [A-z]{3} +[0-9]{1,2} ",
                           "[0-9]{2}:[0-9]{2}:[0-9]{2} [0-9]{4}"),
                    metadata(p)$created))
})

