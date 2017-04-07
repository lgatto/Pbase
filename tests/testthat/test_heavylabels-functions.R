context("heavylabels-functions")

## test_that(".calculateHeavyLabels", {
##     path <- system.file(file.path("extdata"), package="Pbase")

##     p <- Proteins(file.path(path, "heavylabels_proteins_of_interest.fasta"))
##     r <- read.csv(file.path(path, "heavylabels_groundtruth.csv"),
##                 stringsAsFactors = FALSE)

##     ## the fasta file has invalid comments
##     ## names(p@aa) <- ametadata(p)$Comment
##     p@aa <- aa(p)[order(names(aa(p)))]

##     ## cleave the proteins using PLGS rule
##     pc <- cleave(p, custom = "[KR](?=[^P])")

##     peptides <- setNames(r$Peptide, r$Protein)

##     ## calculate labeled peptides
##     heavyLabels <- calculateHeavyLabels(pc, peptides, endsWith = NULL)

##     ## ensure same order
##     i <- match(paste0(r$Protein, r$Peptide),
##                paste0(heavyLabels$Protein, heavyLabels$Peptide))
##     heavyLabels <- heavyLabels[i, ]

##     ## overwrite rownames
##     rownames(heavyLabels) <- NULL

##     ## compare only the following columns
##     cols <- c("Protein", "Peptide", "spikeTideResult", "spikeTide")

##     expect_equal(r[, cols], heavyLabels[, cols])
##     expect_error(calculateHeavyLabels(pc, unname(peptides)),
##                  paste0("No names for ", sQuote("peptides"), " available!"))
##     expect_error(calculateHeavyLabels(p, peptides),
##                  paste0("You have to ", sQuote("cleave"), " your proteins first!"))
## })

