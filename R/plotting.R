plotAsGeneRegionTrack <- function(..., genome = "hg38",
                                  plot = TRUE) {
    args <- pairlist(...)
    stopifnot(sapply(args, is, "GRanges"))
    args <- lapply(args, GeneRegionTrack)
    if (!is.null(genome)) {
        chr <- as.character(seqnames(args[[1]])[1])
        ideoTrack <- IdeogramTrack(genome = genome,
                                   chromosome = chr)
        axisTrack <- GenomeAxisTrack()
        args <- c(ideoTrack, axisTrack, args)
        if (plot) plotTracks(args, add53 = TRUE, add35 = TRUE)
    } else {
        if (plot) plotTracks(args)
    }
    invisible(args)
}


plotAsAnnotationTrack <- function(x, ..., genome = "hg38") {
    stopifnot(is(x, "GRanges"))
    args <- plotAsGeneRegionTrack(..., genome = genome, plot=FALSE)
    pepTrack <- AnnotationTrack(start = start(x),
                                end = end(x),
                                chr = rtracklayer::chrom(x),
                                strand = strand(x),
                                group = mcols(x)$group,
                                id = mcols(x)$pepseq,
                                fill = ifelse(mcols(x)$exonJunctions,
                                    "red", "steelblue"),
                                col = NULL,
                                name = "peptides")
    args <- c(args, pepTrack)
    plotTracks(args, groupAnnotation = "id")
    invisible(args)
}
