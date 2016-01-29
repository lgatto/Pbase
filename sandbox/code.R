library("Pbase")
f0 <- "~/Data2/Thermo_HELA_PRT/swissprot_human_canonical_19_09_12.fasta"
p0 <- Proteins(f0)

p1 <- addIdentificationData(p0, "~/Data2/Thermo_HELA_PRT/Thermo_Hela_PRTC_1.mzid",
                            TRUE)

p2 <- addIdentificationData(p0, "~/Data2/Thermo_HELA_PRT/Thermo_Hela_PRTC_1.mzid",
                            FALSE)

## only 10 proteins, chosen as follows
## lns <- elementNROWS(pranges(p0))
## plns <- elementNROWS(aa(p0))
## i <- which(lns > 20)
## j <- which(lns < 6 & lns > 3)
## c(i[1:5], j[1:5])

f <- "./swissprot_human_canonical_10.fasta"
p <- Proteins(f)
p <- addIdentificationData(p, "~/Data2/Thermo_HELA_PRT/Thermo_Hela_PRTC_1.mzid")

save(p, file = "../data/p.rda")

plot(p)


## at <- ATrack(start = c(250, 480), end = c(320, 520),
##              name = "Annotations")
## pat1 <- ProteinAxisTrack()
## pat2 <- ProteinAxisTrack()
## plotTracks(pat, from = 1, to = 850)

## pat <- ProteinAxisTrack(addNC = TRUE, littleTicks = TRUE)
## seqEx <- aa(p)[[1]]
## plotTracks(trackList = c(pat, st), from = 1, to = 40)


## pat<-ProteinAxisTrack(addNC = TRUE, littleTicks = TRUE)
## data(pep_hxb2)
## hxb2_seq <- metadata(pep_hxb2)$sequence
## st <- ProteinSequenceTrack(sequence = hxb2_seq, name = "env")
## s2 <- paste(strsplit(hxb2_seq, "")[[1]][1:100], collapse = "")
## st2 <- ProteinSequenceTrack(sequence = s2, name = "s2")
## plotTracks(trackList = c(pat, st, st2), from = 1, to = 500)

## png("p1.png")
## plotTracks(trackList = c(pat, st, st2), from = 1, to = 200)
## dev.off()

## png("p2.png")
## plotTracks(trackList = c(pat, st, st2), from = 1, to = 1000)
## dev.off()

library("pRoloc")
## load("~/ResearchProjects/collaborations/fusion-data-analysis/data/pdRes-2014-03-22.rda", verbose = TRUE)
load("~/ResearchProjects/collaborations/fusion-data-analysis/data/svmRes-2014-04-15-top.25.rda", verbose=TRUE)
setStockcol(NULL)
setStockcol(paste0(getStockcol(), "80"))

ptsze <- exp(fData(svmRes.2014.03.27)$svm.scores) - 1


pdf("../Figures/fus.pdf")
plot2D(svmRes.2014.03.27, fcol = "svm")
addLegend(svmRes.2014.03.27, fcol = "svm",
          where = "bottomleft", bty = "n",
          cex = .8)
dev.off()

########################################

library("biomaRt")

## Accession: P53501
## Name: ACT3_DROME
## FBgn0000044

## P60709 (ACTB_HUMAN) 
## ENSG00000075624

ens <- useMart("ensembl", "dmelanogaster_gene_ensembl")

select(ens, keys = "P53501",
       keytype = "uniprot_swissprot",
       columns = c(columns(ens)[1:12],
           "uniprot_swissprot"))

select(ens, keys = "ACT3_DROME",
       keytype = "uniprot_swissprot",
       columns = c(columns(ens)[1:12],
           "uniprot_swissprot"))

ens <- useMart("ensembl", "hsapiens_gene_ensembl")

select(ens, keys = "P60709",
       keytype = "uniprot_swissprot",
       columns = c(columns(ens)[1:12],
           "uniprot_swissprot"))

select(ens, keys = "ACTB_HUMAN",
       keytype = "uniprot_swissprot",
       columns = c(columns(ens)[1:12],
           "uniprot_swissprot"))



ens <- useMart("ensembl", "mmusculus_gene_ensembl")


bm <- select(ens, keys = "DTX4_MOUSE",
             keytype = "uniprot_swissprot",
             columns = c("ensembl_exon_id",
                 "peptide",
                 "exon_chrom_start",
                 "exon_chrom_end"))

bm <- select(ens, keys = "DTX4_MOUSE",
             keytype = "uniprot_swissprot",
             columns = c(
                 "ensembl_gene_id",
                 "ensembl_transcript_id",
                 "ensembl_peptide_id",
                 "ensembl_exon_id",  
                 "description",     
                 "chromosome_name",
                 "strand",
                 "exon_chrom_start",
                 "exon_chrom_end",
                 "5_utr_start",
                 "5_utr_end",
                 "3_utr_start",
                 "3_utr_end",
                 "start_position",     
                 "end_position"))

                 ## "band",
                 ## "transcript_start",
                 ## "transcript_end",
                 ## "exon_chrom_start",
                 ## "exon_chrom_end",
                 ## "rank",
                 ## "uniprot_swissprot",
                 ## "uniprot_swissprot",
                 ## "gene_exon",
                 ## "cdna",
                 ## "coding",
                 ## "peptide"))


library("AnnotationHub")
hub <- AnnotationHub()

## data exploration
length(names(hub))  ## resources available
md <- metadata(hub) ## DataFrame

dna <- hub$ensembl.release.74.fasta.mus_musculus.dna.Mus_musculus.GRCm38.74.dna.toplevel.fa.rz
## dnax <- readDNAStringSet(dna$index)
pep <- hub$ensembl.release.74.fasta.mus_musculus.pep.Mus_musculus.GRCm38.74.pep.all.fa.rz
gtf <- hub$ensembl.release.74.gtf.mus_musculus.Mus_musculus.GRCm38.74.gtf_0.0.1.RData

i <- which(gtf$protein_id == "ENSMUSP00000040229")
gtf[i]

library("Gviz")

library("BSgenome")
library("BSgenome.Mmusculus.UCSC.mm10")
mus <- Mmusculus
## getSeq(mus, gtf[i])

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
xdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
##txdb <- extractTranscriptSeqs(mus, TxDb.Mmusculus.UCSC.mm10.knownGene)


####################
ens <- useMart("ensembl", "mmusculus_gene_ensembl")

## was not mapped by Alejandro
bm <- select(ens, keys = "SARNP_MOUSE",
             keytype = "uniprot_swissprot",
             columns = c(
                 "ensembl_gene_id",
                 "ensembl_transcript_id",
                 "ensembl_peptide_id",
                 "ensembl_exon_id",  
                 "description",     
                 "chromosome_name",
                 "strand",
                 "exon_chrom_start",
                 "exon_chrom_end",
                 "5_utr_start",
                 "5_utr_end",
                 "3_utr_start",
                 "3_utr_end",
                 "start_position",     
                 "end_position"))

###################################


library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

tx_seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19, txdb)

tr <- transcripts(txdb)

(cds <- cdsBy(txdb))
cds_seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19, cds)
