## load Pbase
library("devtools")
library("mzID")
library("cleaver")
devtools::load_all(file.path("..", ".."))

p <- Proteins("proteins_of_interest.fasta")
r <- read.csv("result_26_09_13.csv", stringsAsFactors = FALSE)

## the fasta file has invalid comments
names(p@aa) <- ametadata(p)$Comment
p@aa <- aa(p)[order(names(aa(p)))]

## stupid PLGS rule?
pc <- cleave(p, custom = "[KR](?=[^P])")

peptides <- setNames(r$Peptide, r$Protein)
heavyLabels <- .calculateHeavyLabels(peptides, pc, endsWith = NULL)

## ensure same order
i <- match(paste0(r$Protein, r$Peptide),
           paste0(heavyLabels$Protein, heavyLabels$Peptide))
heavyLabels <- heavyLabels[i, ]

result <- heavyLabels$spikeTide == r$spikeTide

message("Overlap: ", round(mean(result), 2),
        " (", sum(result), "/", length(result), ")")

verified <- c("FVTTDTR", "IDLVIAK", "VLTAYLK", "YFNPAFYLR", "IGGGGEGSWFR")

message("All mismatches verified as correct? : ",
        all(heavyLabels$Peptide[!result] %in% verified))
