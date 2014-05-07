context("fasta-functions")

test_that(".fastaCommentParser", {
  fastacmt <- c(paste0(
    ">tr|B1XC03|B1XC03_ECODH HU, DNA-binding transcriptional regulator, ",
    "alpha subunit OS=Escherichia coli (strain K12 / DH10B) GN=hupA PE=3 SV=1"),
                paste0(
    ">sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha ",
    "chain OS=Homo sapiens GN=HLA-A PE=1 SV=2"),
                paste0(
    ">sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha ",
    "chain OS=Homo sapiens PE=1 SV=2"),
                paste0(
    ">sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha ",
    "OS=Macaca fascicularis GN=YWHAB"),
                paste0(
    ">FOO|NO_ID_HERE|Not a random string."))
  result <- matrix(c(
    "tr", "sp", "sp", "sp", NA_character_,
    "B1XC03", "P04439", "P04439", "Q4R572-2", NA_character_,
    "B1XC03_ECODH", "1A03_HUMAN", "1A03_HUMAN", "1433B_MACFA",
    NA_character_,
    NA_character_, NA_character_, NA_character_, "Short",
    NA_character_,
    "HU, DNA-binding transcriptional regulator, alpha subunit",
    "HLA class I histocompatibility antigen, A-3 alpha chain",
    "HLA class I histocompatibility antigen, A-3 alpha chain",
    "14-3-3 protein beta/alpha", NA_character_,
    "Escherichia coli (strain K12 / DH10B)", "Homo sapiens", "Homo sapiens",
    "Macaca fascicularis", NA_character_,
    "hupA", "HLA-A", NA_character_, "YWHAB", NA_character_,
    "3", "1", "1", NA_character_, NA_character_,
    "1", "2", "2", NA_character_, NA_character_), ncol = 9)
  expect_identical(Pbase:::.fastaCommentParser(fastacmt), result)
})

test_that(".fastaComments2DataFrame", {
  fastacmt <- c(paste0(
    ">tr|B1XC03|B1XC03_ECODH HU, DNA-binding transcriptional regulator, ",
    "alpha subunit OS=Escherichia coli (strain K12 / DH10B) GN=hupA PE=3 SV=1"),
                paste0(
    ">sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha ",
    "chain OS=Homo sapiens GN=HLA-A PE=1 SV=2"),
                paste0(
    ">sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha ",
    "chain OS=Homo sapiens PE=1 SV=2"),
                paste0(
    ">sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha ",
    "OS=Macaca fascicularis GN=YWHAB"),
                paste0(
    ">FOO|NO_ID_HERE|Not a random string"))

  ir <- IRanges(start = 1L:5L, end = 11L:15L)

  df <- DataFrame(
    DB = Rle(factor(c("tr", "sp", "sp", "sp", NA_character_))),
    AccessionNumber = c("B1XC03", "P04439", "P04439", "Q4R572-2", "Pb1"),
    EntryName = c("B1XC03_ECODH", "1A03_HUMAN", "1A03_HUMAN", "1433B_MACFA",
                  NA_character_),
    IsoformName = Rle(c(NA_character_, NA_character_, NA_character_, "Short",
                        NA_character_)),
    ProteinName = c("HU, DNA-binding transcriptional regulator, alpha subunit",
                    "HLA class I histocompatibility antigen, A-3 alpha chain",
                    "HLA class I histocompatibility antigen, A-3 alpha chain",
                    "14-3-3 protein beta/alpha", NA_character_),
    OrganismName = Rle(factor(c("Escherichia coli (strain K12 / DH10B)",
                                "Homo sapiens", "Homo sapiens",
                                "Macaca fascicularis", NA_character_))),
    GeneName = Rle(factor(c("hupA", "HLA-A", NA_character_, "YWHAB",
                            NA_character_))),
    ProteinExistence = Rle(factor(c(3, 1, 1, NA, NA),
                                  labels = c("Evidence at protein level",
                                             "Evidence at transcript level",
                                             "Inferred from homology",
                                             "Predicted",
                                             "Uncertain"),
                                  levels = 1L:5L)),
    SequenceVersion = Rle(c("1", "2", "2", NA_character_, NA_character_)),
    Comment = Rle(c(rep(NA_character_, 4),
                    "FOO|NO_ID_HERE|Not a random string")),
    Filename = Rle(factor(rep("uniprot.fasta", 5L))))

  result <- ir
  mcols(result) <- df

  expect_equal(Pbase:::.addFastaInformation2mcol(ir, fastacmt, "uniprot.fasta"),
               result)
})

test_that(".addFastaInformation2Mcol", {
  fastacmt <- c(paste0(
    ">tr|B1XC03|B1XC03_ECODH HU, DNA-binding transcriptional regulator, ",
    "alpha subunit OS=Escherichia coli (strain K12 / DH10B) GN=hupA PE=3 SV=1"),
                paste0(
    ">sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha ",
    "chain OS=Homo sapiens GN=HLA-A PE=1 SV=2"),
                paste0(
    ">sp|P04439|1A03_HUMAN HLA class I histocompatibility antigen, A-3 alpha ",
    "chain OS=Homo sapiens PE=1 SV=2"),
                paste0(
    ">sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha ",
    "OS=Macaca fascicularis GN=YWHAB"),
                paste0(
    ">FOO|NO_ID_HERE|Not a random string"))

  result = DataFrame(
    DB = Rle(factor(c("tr", "sp", "sp", "sp", NA_character_))),
    AccessionNumber = c("B1XC03", "P04439", "P04439", "Q4R572-2", "Pb1"),
    EntryName = c("B1XC03_ECODH", "1A03_HUMAN", "1A03_HUMAN", "1433B_MACFA",
                  NA_character_),
    IsoformName = Rle(c(NA_character_, NA_character_, NA_character_, "Short",
                        NA_character_)),
    ProteinName = c("HU, DNA-binding transcriptional regulator, alpha subunit",
                    "HLA class I histocompatibility antigen, A-3 alpha chain",
                    "HLA class I histocompatibility antigen, A-3 alpha chain",
                    "14-3-3 protein beta/alpha", NA_character_),
    OrganismName = Rle(factor(c("Escherichia coli (strain K12 / DH10B)",
                                "Homo sapiens", "Homo sapiens",
                                "Macaca fascicularis", NA_character_))),
    GeneName = Rle(factor(c("hupA", "HLA-A", NA_character_, "YWHAB",
                            NA_character_))),
    ProteinExistence = Rle(factor(c(3, 1, 1, NA, NA),
                                  labels = c("Evidence at protein level",
                                             "Evidence at transcript level",
                                             "Inferred from homology",
                                             "Predicted",
                                             "Uncertain"),
                                  levels = 1L:5L)),
    SequenceVersion = Rle(c("1", "2", "2", NA_character_, NA_character_)),
    Comment = Rle(c(rep(NA_character_, 4),
                    "FOO|NO_ID_HERE|Not a random string")))
  expect_equal(Pbase:::.fastaComments2DataFrame(fastacmt), result)

})

test_that(".isUniProtAccessionNumber", {
  acn <- c("B1XC03", "P04439", "P04439", "Q4R572-2", ## valid
           "Pb123", "A12ABC0", "A1za111", "P00-34")  ## not valid
  result <- rep(c(TRUE, FALSE), each = 4L)
  expect_identical(Pbase:::.isUniProtAccessionNumber(acn), result)
})

test_that(".isPbaseAccessionNumber", {
  acn <- c("Pb1", "Pb1000", "PB10", "B1XC03")
  result <- rep(c(TRUE, FALSE), each = 2L)
  expect_identical(Pbase:::.isPbaseAccessionNumber(acn), result)
})

test_that(".isValidAccesionNumber", {
  acn <- c("B1XC03", "P04439", "P04439", "Q4R572-2", ## valid
           "Pb123", "A12ABC0", "A1za111", "P00-34")  ## not valid
  result <- c(rep(TRUE, 5L) , rep(FALSE, 3))
  expect_identical(Pbase:::.isValidAccessionNumber(acn), result)
})

