## ----style, echo = FALSE, results = 'asis', message=FALSE----------------
BiocStyle::markdown()

## ----bm------------------------------------------------------------------
tr <- "ENST00000373316"

## ----, ens---------------------------------------------------------------
suppressMessages(library("Gviz"))
suppressMessages(library("biomaRt"))

ens <- useMart("ensembl", "hsapiens_gene_ensembl")
etr <- BiomartGeneRegionTrack(biomart = ens,
                              transcript = tr,
                              genome = "hg19")
etr <- split(etr, transcript(etr))
etr <- etr[[tr]]
etr <- ranges(etr)

## ----ucsc----------------------------------------------------------------
utr <- BiomartGeneRegionTrack(transcript = tr,
                              genome = "hg19")
utr <- split(utr, transcript(utr))
utr <- utr[[tr]]
utr <- ranges(utr)

## ----chain---------------------------------------------------------------
library("AnnotationHub")
hub <- AnnotationHub()
query(hub, 'hg19ToHg38')
chain <- query(hub, 'hg19ToHg38')[[1]]

## ------------------------------------------------------------------------
library("rtracklayer")

res <- liftOver(utr, chain)
res <- unlist(res)

## set annotation
genome(res) <- "hg19"
names(res) <- NULL

identical(res, etr)

## ----si------------------------------------------------------------------
sessionInfo()

