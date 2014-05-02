library("Pbase")

fas <- "swissprot_human_canonical_19_09_12.fasta"

p <- Proteins(fas)[1:10]
pfeatures(p)
pp <- cleave(p[1:10])
pfeatures(pp)
