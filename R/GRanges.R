##' This function takes on or more Ensembl transcript identifiers,
##' queries Biomart and constructs a \code{GRangesList} object as
##' would \code{Gviz::BiomartGeneRegionTrack} for a genomic region (in
##' fact, currently most of the code has been taken from
##' \code{Gviz::.fetchBMData} and \code{GViz::.chrName} is used to
##' validate chromosome names).
##'
##' @title From a transcript identifier to \code{GRanges} object
##' @param etrid A vector of Ensembl transcript identifiers.
##' @param ens A instance of class \code{Mart} from biomaRt. If
##' missing, \code{useMart("ensembl", "hsapiens_gene_ensembl")} is
##' used.
##' @param use.names If set to \code{TRUE} and \code{etrid} has names,
##' then the latter are used to name the output.
##' @return A \code{GRangesList} object of length
##' \code{length(etrid)}.
##' @author Laurent Gatto
##' @examples
##' id <- c("ENST00000612959", "ENST00000317091")
##' grl1 <- etrid2grl(id[1])
##' grl1
##' grl <- etrid2grl(id)
##' stopifnot(all.equal(id, names(grl)))
etrid2grl <- function(etrid, ens, use.names = FALSE) {    
    if (missing(ens))
        ens <- useMart("ensembl", "hsapiens_gene_ensembl")
    if (!validObject(ens))
        stop("The Mart instance is not valid.")
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

    gr <- GRanges(seqnames=bm$chromosome,
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
    gr <- sort(gr)
    grl <- split(gr, mcols(gr)$transcript)
    grl <- grl[etrid] ## same order as input ids

    if (use.names & !is.null(names(etrid)))
        names(grl) <- names(etrid)

    if (validObject(grl))
        return(grl)
    
}


setMethod("proteinCoding", "GRanges",
          function(object, mcol = "feature",
                   coding = "protein_coding") {
              sel <- mcols(object)[, mcol] == coding
              object[sel, ]
          })

setMethod("proteinCoding", "GRangesList",
          function(object, mcol = "feature", coding = "protein_coding")
              endoapply(object, proteinCoding))

