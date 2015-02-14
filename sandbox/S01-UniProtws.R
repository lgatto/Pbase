library("UniProt.ws")

head(availableUniprotSpecies())

availableUniprotSpecies(pattern="sapiens")
#   taxon ID                  Species name
# 1    63221 Homo sapiens neanderthalensis
# 2     9606                  Homo sapiens

## Homo sapiens is the default
taxId(UniProt.ws) <- 9606

## available keys
keytypes(UniProt.ws)

## available columns
columns(UniProt.ws)

## ALBU_HUMAN
keys <- c("P02768")
keytype <- "UNIPROTKB"
cols <- c("UNIPROTKB", "ID", "LENGTH",  "SEQUENCE")

res <- select(UniProt.ws, keys, cols, keytype)
