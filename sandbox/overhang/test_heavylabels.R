## load Pbase
library("devtools")
devtools::load_all(file.path("..", ".."))

p <- Proteins("proteins_of_interest.fasta")
r <- read.csv("result_26_09_13.csv", stringsAsFactors = FALSE)

## the fasta file has invalid comments
names(p@aa) <- ametadata(p)$Comment

## stupid PLGS rule?
pc <- cleave(p, custom="[KR](?=[^P])")

cols <- c("Peptide", "N_overhang", "C_overhang", "spikeTideResult", "spikeTide")

result <- logical(nrow(r))

for (id in seqnames(pc)) {
    for (i in which(r$Protein == id)) {
        d <- dummyAddOverhang(setNames(r$Peptide[i], id), pc, endsWith = NULL)

        spikeTide <- paste0(na.omit(c(d$N_overhang, d$Peptide, d$C_overhang)),
                            collapse=".")

        if (r$spikeTideResult[i] == d$spikeTideResult &&
            r$spikeTide[i] == spikeTide) {
            result[i] <- TRUE
        } else {
            message("\n ", r$Peptide[i], " NOT equal")
            message("Pavel's list:")
            print(r[i, cols])
            message("Pbase result:")
            print(d)
            print(spikeTide)
        }
    }
}

message("Overlap: ", round(mean(result), 2),
        " (", sum(result), "/", length(result), ")")
