library("Biostrings")
library("cleaver")

setClass("ProteinsIRL",
         slots = list(
          aa = "AAStringSet",
          pfeatures = "CompressedIRangesList"))

setClass("ProteinsAASSL",
         slots = list(
          aa = "AAStringSet",
          pfeatures = "AAStringSetList"))

pir <- new("ProteinsIRL")
pal <- new("ProteinsAASSL")

f <- readAAStringSet("swissprot_human_canonical_19_09_12.fasta")

tracemem(f)
pir@aa <- f
pal@aa <- f

## no copy
tracemem(pir@aa)
tracemem(pal@aa)

system.time(
pir@pfeatures <- cleavageRanges(pir@aa)
)
system.time(
pal@pfeatures <- cleave(pal@aa)
)

object.size(pir)
object.size(pal)

str(pir)
str(pal)

## here the data are copied
system.time(
pirpetides <- extractAt(pir@aa, unname(pir@pfeatures))
)

## here not
system.time(
palpeptides <- pal@pfeatures
)
