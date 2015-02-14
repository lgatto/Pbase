### TODO:
### - what happens when there is no mapping
### - consider a mapToGenome,Proteins,GRangesList where
###   all proteins are mapped against all ranges.
### - unit tests

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
            accession = seqnames(pObj)[1],
            gene = mcols(grObj)$gene[1],
            transcript = mcols(grObj)$transcript[1],
            symbol = mcols(grObj)$symbol[1],
            exonJunctions = junc)    
}    


setMethod("mapToGenome", c("Proteins", "GenomicRanges"),
          function(x, genome, ...) .mapToGenome(x, genome, ...))


setMethod("pmapToGenome", c("Proteins", "GRangesList"),
          function(x, genome, ...) {
              if (length(x) != length(genome))
                  stop("'x' and 'genome' must have the same length")
              
              l <- vector("list", length = length(x))
              names(l) <- names(genome)
              for (i in seq_len(length(x))) 
                  l[[i]] <- .mapToGenome(x[i], genome[[i]])
              ans <- GRangesList(l)
              if (validObject(ans))
                  return(ans)
          })

