##' These functions convert ranges of peptides or exons to
##' \code{AnnotationTrack} or \code{GeneRegionTrack} objects from the
##' \code{Gviz} package and produces the corresponding plot. The
##' \code{genome} argument controls whether additional ideogram and
##' axis tracks are to be plotted. \code{plotAsAnnotationTrack} plots
##' peptides that span multiple exons in red and connects them with a
##' grey line.
##'
##' @title Plot gene region and annotation tracks
##' @param x A \code{Granges} object containing peptides genomics
##'     coordinates. These ranges are converted to a \code{AnnotationTrack}.
##' @param ... One or more \code{GRanges} instances, typically resulting from
##'     calling \code{\link{etrid2grl}}, or, a single \code{GRangesList}. These
##'     ranges are converted to \code{GeneRegionTrack} instances.
##' @param genome A \code{character} of length 1, giving the name of the
##'     genome. Default is \code{"hg38"}. If \code{NULL}, no chromosome and axis
##'     tracks are displayed.
##' @param plot A \code{logical} defining if the figure should be
##'     plotted. Default is \code{TRUE}.
##' @return Used for its plotting side effects. Invisible returns a list of
##'     tracks.
##' @author Laurent Gatto
##' @aliases plotAsGeneRegionTrack
plotAsAnnotationTrack <- function(x, ..., genome = "hg38",
                                  plot = TRUE) {
    stopifnot(is(x, "GRanges"))
    args <- plotAsGeneRegionTrack(..., genome = genome, plot=FALSE)
    pepTrack <- AnnotationTrack(x,
                                group = mcols(x)$group,
                                id = mcols(x)$pepseq,
                                col = NULL,
                                groupAnnotation = "id",
                                name = "peptides")
    feature(pepTrack) <- ifelse(mcols(x)$exonJunctions,
                                "exex", "ex")
    args <- c(args, pepTrack)
    if (plot) plotTracks(args, exex = "red", ex = "steelblue")
    invisible(args)
}

##' @rdname plotAsAnnotationTrack
plotAsGeneRegionTrack <- function(..., genome = "hg38",
                                  plot = TRUE) {
    args <- pairlist(...)
    if (length(args) == 1 & is(args[[1]], "GRangesList"))
        args <- args[[1]]
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
