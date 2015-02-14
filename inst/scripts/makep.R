gzmgf <- "../extdata/Thermo_Hela_PRTC_selected.mgf.gz"
mgf <- R.utils::gunzip(gzmgf)
fas <- "../extdata/HUMAN_2015_02_selected.fasta"
msgfpath <- "~/bin/MSGFPlus.20140630/MSGFPlus.jar"
mzid <- "Thermo_Hela_PRTC_selected.mzid"

msgfrun <- paste("java -Xmx3500M -jar", msgfpath,
                 "-s", mgf, "-d", fas, "-o", mzid)

system(msgfrun)

library("Pbase")
p <- Proteins(fas)
p <- addIdentificationData(p, mzid)

## We also want ENST for the mapping vignette.  Adding manually
## here.

## ENST 1 retrieved manually from Ensembl
## ENST 2:9 identified through the
## README_HUMAN_uniprot_ensembl_gene_coordindate.txt coordinate file.

ENST <- c(A4UGR9 = "ENST00000409195",
          A6H8Y1 = "ENST00000358731",
          O43707 = "ENST00000252699",
          O75369 = "ENST00000295956",
          P00558 = "ENST00000373316",
          P02545 = "ENST00000368300",
          P04075 = "ENST00000338110",
          'P04075-2' = "ENST00000395248",
          P60709 = "ENST00000331789")
stopifnot(identical(seqnames(p), names(ENST)))
mcols(p@aa)$ENST <- ENST
metadata(p, "Species") <- "Homo sapiens"
metadata(p, "UniProtRelease") <- "2015_02"
metadata(p, "Genome") <- "GRCh38"

stopifnot(validObject(p))

## Saving to extdata and data
file.rename(mzid,
            file.path("../extdata", mzid))

save(p, file = "../../data/p.rda",
     compress = "xz", compression_level = 9)

library("MSnbase")
pms <- readMgfData(mgf)
stopifnot(validObject(pms))
save(pms, file = "../../data/pms.rda",
     compress = "xz", compression_level = 9)

## gzip mgf file
R.utils::gzip(mgf)

## cleaning up

unlink(file.path("../extdata",
                 c("HUMAN_2015_02_selected.canno",
                   "HUMAN_2015_02_selected.cnlcp",
                   "HUMAN_2015_02_selected.csarr",
                   "HUMAN_2015_02_selected.cseq")))
