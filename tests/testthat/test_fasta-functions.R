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
    "OS=Macaca fascicularis GN=YWHAB"))
  result = matrix(c(
    "tr", "sp", "sp", "sp",
    "B1XC03", "P04439", "P04439", "Q4R572-2",
    "B1XC03_ECODH", "1A03_HUMAN", "1A03_HUMAN", "1433B_MACFA",
    NA_character_, NA_character_, NA_character_, "Short",
    "HU, DNA-binding transcriptional regulator, alpha subunit",
    "HLA class I histocompatibility antigen, A-3 alpha chain",
    "HLA class I histocompatibility antigen, A-3 alpha chain",
    "14-3-3 protein beta/alpha",
    "Escherichia coli (strain K12 / DH10B)", "Homo sapiens", "Homo sapiens",
    "Macaca fascicularis",
    "hupA", "HLA-A", NA_character_, "YWHAB",
    "3", "1", "1", NA_character_,
    "1", "2", "2", NA_character_), ncol = 9)
  expect_identical(Pbase:::.fastaCommentParser(fastacmt), result)
})

