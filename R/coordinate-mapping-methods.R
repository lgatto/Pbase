setGeneric("mapToGenome",
           function(x, genome, ...) standardGenerics("mapToGenome"))

setGeneric("pmapToGenome",
           function(x, genome, ...) standardGenerics("pmapToGenome"))


##' TODO
##' 
##' To calculate the genomic coordinates of peptides we
##' 
##' 1) Calculate the contiguous exon coordinates on the cDNA sequence.
##' 2) Calculate the peptide coordinates on the cDNA sequence from the
##'    original peptide ranges along the protein sequence.
##' 3) We calculate the indices of the exons in which all the peptides
##'    start and end on the cDNA sequence.
##' 4) We infer the coordinates on the genome from the actual peptide
##'    starting/ending position and exon index on the cDNA sequence, the
##'    exons coordinates on the cDNA sequences and the matching exons
##'    coordinates on the genome.
##'
##' @title Map peptides on genomic coordinates
##' @param pObj An \code{Proteins} instance of length 1, contain
##' \code{pfeatures} and \code{pranges}.
##' @param grObj A \code{GRanges} object containing the exon positions
##' of the protein in \code{pObj}, as created by
##' \code{\link{etrid2grl}}.
##' @param ... Other argumants
##' @return A \code{Granges} object with the peptide coordinates on
##' the genome.
##' @seealso The \code{\link{proteinCoding}} function to remove
##' non-protein coding ranges before maping peptides to their genomic
##' coordinates.
##' @noR
.mapToGenome <- function(pObj, grObj, ...) {
    ## exons ranges along the protein (1)
    j <- cumsum(width(grObj))
    i <- cumsum(c(1, width(grObj)))[1:length(j)]
    prex <- IRanges(start = i, end = j)

    ## peptide position on protein
    peprngProt <- pranges(pObj)[[1]]
    ## peptide positions on cdna (2)
    peprngCdna <- IRanges(start = 1 + (start(peprngProt)-1) * 3,
                          width = width(peprngProt) * 3)


    ## find exon and position in prex (3)
    start_ex <- subjectHits(findOverlaps(start(peprngCdna), prex))
    end_ex <- subjectHits(findOverlaps(end(peprngCdna), prex))

    if (any(junc <- start_ex != end_ex))
        message("Peptide(s) ", paste(which(junc), collapse = ", "),
                " overlap(s) exon junctions.")    

    ## (4)
    getPos <- function(pos, idx, nclex, prtex) {
        ## position in cdna
        ## exon index
        start(nclex[idx]) + (pos - start(prtex[idx]))
    }
    peptides_on_genome <-
        IRanges(start = getPos(start(peprngCdna), start_ex, grObj, prex),
                end = getPos(end(peprngCdna), end_ex, grObj, prex))

    chr <- as.character(seqnames(grObj)@values)
    
    GRanges(seqnames = rep(chr, length(peptides_on_genome)),
            ranges = peptides_on_genome,
            strand = strand(grObj)@values,
            pepseq = as.character(pfeatures(pObj)[[1]]),
            accession = seqnames(pObj),
            gene = mcols(grObj)$gene[1],
            transcript = mcols(grObj)$transcript[1],
            symbol = mcols(grObj)$symbol[1],
            exonJunctions = junc)    
}    

