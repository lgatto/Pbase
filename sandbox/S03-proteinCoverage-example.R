library("Pbase")
library("cleaver")
library("Biostrings")

p <- Proteins("../inst/extdata/04_test_database.fasta")

aa <- cleave(aa(p))[[1]]
mcols(aa) <- ametadata(p)[rep(1, 14),]

set.seed(1)
aa <- aa[sample(length(aa), size = 5)]

## it is horrible slow!
Rprof()
system.time(
pvs1 <- proteinCoverage(p, aa)
)
Rprof(NULL)
summaryRprof()

pvs1

Pbase:::.coverageOutputTextAAStringSet(pvs1@aa, pvs1@coverage)[1]
