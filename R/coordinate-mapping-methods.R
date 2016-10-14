### gr: and GRanges object with mapped peptides
### j: a logical, typically mcols(gr)$exonJunctions
### ex: a GRanges object with exonsx
splitExonJunctions <- function(gr, j, ex) {
    ## (1) the ranges to be split
    gr2split <- gr[j]

    ## (2) remove junc exons from
    gr <- gr[!j]

    ## (3) split ranges, but missing mcols
    ## this must be done range by range, as we want to preserve
    ## duplicated ranges and ranges that have the same start in exon i
    ## and different ends in exon i+1, and would be split as 3 splits.
### THIS IS VERY SLOW
    grsplit <- sapply(gr2split, intersect, ex)
    for (i in 1:length(grsplit))
         mcols(grsplit[[i]])$._N_ <- i
    grsplit <- Reduce(c, grsplit)

    ## (4) add mcols
    groups <- mcols(grsplit)$._N_
    mcols(grsplit) <- mcols(gr2split)[mcols(grsplit)$._N_, ]
    mcols(grsplit)$group <- groups

    if (length(gr) > 0)
        mcols(gr)$group <- (max(groups)+1):(max(groups)+length(gr))

    ## (5) add back to original ranges
    ans <- c(gr, grsplit)

    sort(ans)
}

setGeneric("mapToGenome",
           function(x, genome, ...) standardGeneric("mapToGenome"))

setGeneric("pmapToGenome",
           function(x, genome, ...) standardGeneric("pmapToGenome"))


## To calculate the genomic coordinates of peptides we
##
## 1) Calculate the contiguous exon coordinates on the cDNA sequence.
## 2) Calculate the peptide coordinates on the cDNA sequence from the
##    original peptide ranges along the protein sequence.
## 3) We calculate the indices of the exons in which all the peptides
##    start and end on the cDNA sequence.
## 4) We infer the coordinates on the genome from the actual peptide
##    starting/ending position and exon index on the cDNA sequence, the
##    exons coordinates on the cDNA sequences and the matching exons
##    coordinates on the genome.
##

tryCatchMapToGenome <- function(pObj, grObj, ...)
    tryCatch(.mapToGenome(pObj, grObj, ...),
             warning = function() NULL,
             error = function(e) {
                 warning("Mapping failed. Returning an empty range.",
                         call. = FALSE)
                 GRanges()
             })

### 'pObj' is a Proteins of length 1
### 'grObj' is a GRanges
### Returns a GRanges with positions of pranges(pObj)[[1]] mapped to grObj
.mapToGenome <- function(pObj, grObj, ...) {
    if (length(pObj) > 1) {
        warning("Only considering first protein in the Proteins object.")
        pObj <- pObj[1]
    }

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
    ## In which exon does the i-th peptide start?
    start_ex <- subjectHits(findOverlaps(start(peprngCdna), prex))
    ## In which exon does the i-th peptide end?
    end_ex <- subjectHits(findOverlaps(end(peprngCdna), prex))

    junc <- start_ex != end_ex

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

    x <- GRanges(seqnames = rep(chr, length(peptides_on_genome)),
                 ranges = peptides_on_genome,
                 strand = strand(grObj)@values,
                 pepseq = as.character(pfeatures(pObj)[[1]]),
                 accession = seqnames(pObj)[1],
                 gene = mcols(grObj)$gene[1],
                 transcript = mcols(grObj)$transcript[1],
                 symbol = mcols(grObj)$symbol[1],
                 exonJunctions = junc)

    if (any(mcols(x)$exonJunctions))
        x <- splitExonJunctions(x, mcols(x)$exonJunctions, grObj)
    else
        mcols(x)$group <- 1:length(x)
    if (validObject(x))
        return(x)
}

############################################################
## .mapToGenome2
##' maps ranges within the protein sequence to genomic sequence.
##' This function ensures:
##' 1) that the provided grObj could indeed encode the protein in pObj (based on
##'    its sequence length).
##' 2) that coordinates are correctly returned for + and - strand genes.
##' @param pObj is a Proteins of length 1
##' @param grObj is a GRanges representing the start and end coordinates of the
##' (CDS!) exons of the transcript encoding the protein sequence.
##' @return Returns a GRanges with positions of pranges(pObj)[[1]] mapped to
##' grObj
##' @noRd
.mapToGenome2 <- function(pObj, grObj, ...) {
    if (length(pObj) > 1) {
        warning("Only considering first protein in the Proteins object.")
        pObj <- pObj[1]
    }
    ## Check that the length of the CDS corresponds to the length of the AA:
    ## AA length should be (length(CDS) - 1) / 3, -1 because stop codon is not
    ## encoded. Exception: if there is no 3'UTR the length can also be
    ## length(CDS) / 3. For most CDS/AA sequences it works, but there are some,
    ## e.g. ENST00000371584 that have a truncated 5' sequences (no START codon).
    p_width <- width(pObj@aa)
    cds_width <- sum(width(grObj))
    if (p_width != cds_width/3)
        stop("The protein sequence length of the provided Proteins object does",
             " not match the length of the coding sequence in the GRanges!")


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
    ## In which exon does the i-th peptide start?
    start_ex <- subjectHits(findOverlaps(start(peprngCdna), prex))
    ## In which exon does the i-th peptide end?
    end_ex <- subjectHits(findOverlaps(end(peprngCdna), prex))

    junc <- start_ex != end_ex

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

    x <- GRanges(seqnames = rep(chr, length(peptides_on_genome)),
                 ranges = peptides_on_genome,
                 strand = strand(grObj)@values,
                 pepseq = as.character(pfeatures(pObj)[[1]]),
                 accession = seqnames(pObj)[1],
                 gene = mcols(grObj)$gene[1],
                 transcript = mcols(grObj)$transcript[1],
                 symbol = mcols(grObj)$symbol[1],
                 exonJunctions = junc)

    if (any(mcols(x)$exonJunctions))
        x <- splitExonJunctions(x, mcols(x)$exonJunctions, grObj)
    else
        mcols(x)$group <- 1:length(x)
    if (validObject(x))
        return(x)
}


## setMethod("mapToGenome", c("Proteins", "GenomicRanges"),
##           function(x, genome, ...) .mapToGenome(x, genome, ...))

setMethod("pmapToGenome", c("Proteins", "GRangesList"),
          function(x, genome, drop.empty.ranges = TRUE, ...) {
              if (length(x) != length(genome))
                  stop("'x' and 'genome' must have the same length")

              l <- vector("list", length = length(x))
              names(l) <- names(genome)
              for (i in seq_len(length(x)))
                  l[[i]] <- tryCatchMapToGenome(x[i], genome[[i]], ...)
              ans <- GRangesList(l)
              if (drop.empty.ranges)
                  ans <- ans[elementNROWS(ans) > 0]
              if (validObject(ans))
                  return(ans)
          })


setMethod("mapToGenome", c("Proteins", "GRangesList"),
          function(x, genome, drop.empty.ranges = TRUE, ...) {
              if (length(x) == 1 & length(genome) == 1) {
                  ans <- tryCatchMapToGenome(x, genome[[1]], ...)
                  if (drop.empty.ranges & length(ans) == 1)
                      ans <- GRangesList()
                  else {
                      ans <- GRangesList(ans)
                      names(ans) <- names(genome)
                  }
              } else {
                  ## Proteins[n] and genome[m] and mapp all against all
                  ## WITH matching names.

                  nmsx0 <- seqnames(x)
                  nmsg0 <- names(genome)

                  if (is.null(nmsx0) | is.null(nmsg0))
                      stop("'x' and 'genome' must have names.")

                  nmsi <- intersect(nmsx0, nmsg0)
                  ## update input with common names
                  x <- x[nmsx0 %in% nmsi]
                  nmsx <- seqnames(x)
                  genome <- genome[nmsg0 %in% nmsi]

                  if (!all(nmsx0 %in%nmsi))
                      message("Mapping ", length(x), " out of ",
                              length(nmsx0), " peptide ranges.")

                  k <- match(names(genome), seqnames(x))
                  x <- x[k]
                  ans <- pmapToGenome(x, genome, drop.empty.ranges, ...)
                  ## get original seqnames order
                  ans <- ans[order(match(names(ans), nmsx))]
              }

              if (validObject(ans))
                  return(ans)

          })
