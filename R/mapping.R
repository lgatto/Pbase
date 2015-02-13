##' This function takes on or more Ensembl transcript identifiers,
##' queries Biomart and constructs a \code{GRanges} object as would
##' \code{Gviz::BiomartGeneRegionTrack} for a genomic region (in fact,
##' currently most of the code has been taken from
##' \code{Gviz::.fetchBMData} and \code{GViz::.chrName} is used to
##' validate chromosome names).
##' 
##'
##' @title From a transcript identifier to \code{GRanges} object
##' @param etrid A vector of Ensembl transcript identifiers.
##' @param ens A instance of class \code{Mart} from biomaRt. If
##' missing, \code{useMart("ensembl", "hsapiens_gene_ensembl")} is
##' used.
##' @return A \code{Granges} object.
##' @author Laurent Gatto
##' @examples
##' id <- c("ENST00000612959", "ENST00000317091")
##' etrid2gr(id[1])
##' (gr <- etrid2gr(id))
##' ## make it a GRangesList
##' (grl <- split(gr, mcols(gr)$transcript))
etrid2gr <- function(etrid, ens) {    
    if (missing(ens))
        ens <- useMart("ensembl", "hsapiens_gene_ensembl")
    
    bm <- select(ens, keys = etrid,
                 keytype = "ensembl_transcript_id",
                 columns = c(
                     "ensembl_gene_id", "ensembl_transcript_id",
                     "ensembl_exon_id", "exon_chrom_start",
                     "exon_chrom_end", "rank",
                     "strand", "external_gene_name", "gene_biotype",
                     "chromosome_name", "5_utr_start", "5_utr_end",
                     "3_utr_start", "3_utr_end", "phase"))

    ## From GViz:::.fetchBMData
    hasUtr <- !is.na(bm$'5_utr_start') | !is.na(bm$'3_utr_start')
    bmUtr <- bm[hasUtr,, drop=FALSE]
    bmUtr$ffeature <- ifelse(is.na(bmUtr$'5_utr_start'), "utr3", "utr5")
    bmUtr$us <- ifelse(bmUtr$ffeature=="utr3", bmUtr$'3_utr_start', bmUtr$'5_utr_start')
    bmUtr$ue <- ifelse(bmUtr$ffeature=="utr3", bmUtr$'3_utr_end', bmUtr$'5_utr_end')
    bmUtr$'5_utr_end' <- bmUtr$'5_utr_start' <-
        bmUtr$'3_utr_end' <- bmUtr$'3_utr_start' <- NULL
    allUtr <- bmUtr$us == bmUtr$'exon_chrom_start' & bmUtr$ue == bmUtr$'exon_chrom_end'
    utrFinal <- bmUtr[allUtr,, drop=FALSE]
    bmUtr <- bmUtr[!allUtr,, drop=FALSE]
    bmUtrS <- split(bmUtr, ifelse(bmUtr$'exon_chrom_start' == bmUtr$us, "left", "right"))

    utrFinal <- rbind(utrFinal, do.call(rbind, lapply(names(bmUtrS), function(i){
        y <- bmUtrS[[i]]
        if(nrow(y)==0)
            return(NULL)
        yy <- y[rep(1:nrow(y), each=2),]
        sel <- seq(1, nrow(yy), by=2)
        yy[sel, "exon_chrom_end"] <- if(i=="left") yy[sel, "ue"] else yy[sel, "us"]-1
        yy[sel, "ffeature"] <-  yy[sel, ifelse(i=="left", "ffeature", "gene_biotype")]
        yy[sel, "phase"] <-  if(i=="left") -1 else 0
        sel <- seq(2, nrow(yy), by=2)
        yy[sel, "exon_chrom_start"] <- if(i=="left") yy[sel, "ue"]+1 else yy[sel, "us"]
        yy[sel, "ffeature"] <-  yy[sel, ifelse(i=="left", "gene_biotype", "ffeature")]
        yy[sel, "phase"] <- if(i=="left") yy[sel, "phase"] else -1
        yy
    })))

    utrFinal$feature <- utrFinal$ffeature
    bm$feature <- bm$gene_biotype

    keep <- c("ensembl_gene_id","ensembl_transcript_id",
              "ensembl_exon_id","exon_chrom_start",
              "exon_chrom_end", "rank", "strand", "external_gene_name",
              "feature", "chromosome_name", "phase")
    bm <- rbind(bm[!hasUtr,keep, drop=FALSE], utrFinal[,keep])
    bm$chromosome <- Gviz:::.chrName(bm$chromosome, force=TRUE)

    range <- GRanges(seqnames=bm$chromosome,
                     ranges=IRanges(
                         start=bm$'exon_chrom_start',
                         end=bm$'exon_chrom_end'),
                     strand=bm$strand,
                     feature=as.character(bm$feature),
                     gene=as.character(bm$'ensembl_gene_id'),
                     exon=as.character(bm$'ensembl_exon_id'),
                     transcript=as.character(bm$'ensembl_transcript_id'),
                     symbol=as.character(bm$'external_gene_name'),
                     rank=as.numeric(bm$rank),
                     phase=as.integer(bm$phase))
    sort(range)
}
