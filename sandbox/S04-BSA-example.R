## see also: http://proteomicsnews.blogspot.co.uk/
library("Pbase")
library("IRanges")

download.file("http://www.uniprot.org/uniprot/P02769.fasta", "P02769.fasta")
bsa <- Proteins("P02769.fasta")
bsa <- cleave(bsa)
bsa <- pfilter(bsa, mass = c(400, Inf))

## filter signal peptides
bsa@pranges@unlistData@elementMetadata$Filtered[1:4] <- TRUE

## create filtered pranges
bsa@pranges <- IRangesList(P02769=bsa@pranges[[1]][!pcols(bsa)[[1]]$Filtered])

acols(proteinCoverage(bsa))$Coverage
