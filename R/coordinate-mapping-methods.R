### gr: and GRanges object with mapped peptides
### j: a logical, typically mcols(gr)$exonJunctions
### ex: a GRanges object with exonsx
## Update jo: ensure that names are not dropped for spliced features.
## 2nd Update jo: the .mapToGenome2 function does no longer need this function
## as it uses an eventually faster implementation without sapply and for.
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
    mcols(grsplit) <- mcols(gr2split)[groups, ]
    mcols(grsplit)$group <- groups
    ## (4.2) add names
    if (!is.null(names(gr2split)))
        names(grsplit) <- names(gr2split)[groups]

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

tryCatchMapToGenome <- function(pObj, grObj, pcol, ...)
    tryCatch(.mapToGenome2(pObj, grObj, pcol, ...),
             warning = function() NULL,
             error = function(e) {
                 warning("Mapping failed. Returning an empty range.",
                         " Last message was: ", e,
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
##' maps ranges within the protein sequence to genomic sequence. The function
##' takes the acols from the pObj and adds it to the mcols of the returned
##' GRanges. In contrast to the "old" function this:
##' a) does not require the presence of certain columns in the grObj, but uses
##'    all of the mcols that are associated to genes or transcripts.
##' b) Ensure that coordinate conversion for - strand transcripts works.
##' c) Doesn't need splitExonJunctions and thus might be faster.
##' d) Does not throw an error but returns GRanges() if no peptide features
##'    available to map.
##' @param pObj is a Proteins of length 1
##' @param grObj is a GRanges representing the start and end coordinates of the
##' (CDS!) exons of the transcript encoding the protein sequence. Exons HAVE to
##' be ordered by exon index, i.e. first exon is first element, last exon last
##' one. For + strand, exons should be ordered by start increasingly, for -
##' strand by end decreasingly (or start * strand increasingly).
##' @param pcol character(1) defining the column of the \code{pranges}
##' \code{DataFrame} in which the features to be aligned can be found. If
##' missing \code{pvarLabels(pObj)[1]} will be used.
##' @return Returns a GRanges with positions of pranges(pObj)[[1]] mapped to
##' grObj
##' @noRd
.mapToGenome2 <- function(pObj, grObj, pcol, ...) {
    if (length(pObj) > 1) {
        warning("Only considering first protein in the Proteins object.")
        pObj <- pObj[1]
    }
    ## Check Proteins object var presence of features.
    if (length(pcols(pObj)) == 0)
        stop("'pcols' of the provided Proteins object is empty!")
    if (missing(pcol) || length(pcol) == 0) {
        pcol <- pvarLabels(pObj)[1]
    }
    if (!(pcol %in% pvarLabels(pObj)))
        stop("No column named '", pcol, "' present in 'pcol'")
    ## DROP THAT TEST FOR NOW - eventually implement later.
    ## ## Check that the length of the CDS corresponds to the length of the AA:
    ## ## AA length should be (length(CDS) - 1) / 3, -1 because stop codon is not
    ## ## encoded. Exception: if there is no 3'UTR the length can also be
    ## ## length(CDS) / 3. For most CDS/AA sequences it works, but there are some,
    ## ## e.g. ENST00000371584 that have a truncated 5' sequences (no START codon).
    ## ## TODO @jo: fix this; eventually just show a warning.
    ## p_width <- width(pObj@aa)
    ## cds_width <- sum(width(grObj))
    ## if (!(p_width == cds_width/3 | p_width == (cds_width/3 - 1)))
    ##     warning("The protein sequence length of the provided Proteins object",
    ##             " does not match the length of the coding sequence in the",
    ##             " GRanges!")

    strand_num <- 1
    if (as.character(strand(grObj)[1]) == "-") {
        strand_num <- -1
    }
    ## Ensure exons in grObj are ordered correctly.
    exn_order <- order(start(grObj) * strand_num)
    if (!all(exn_order == 1:length(exn_order)))
        stop("Provided exons in 'grObj' are not ordered by exon index/rank!")

    ## exons ranges along the protein (1)
    j <- cumsum(width(grObj))
    i <- cumsum(c(1, width(grObj)))[1:length(j)]
    prex <- IRanges(start = i, end = j)

    ## peptide position on protein
    ## peprngProt <- pranges(pObj)[[1]]
    peprngProt <- pranges(pObj)[[pcol]][[1]]
    if (length(peprngProt) == 0)
        return(GRanges())
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
    ## Translate position within cDNA to position within the respective exon.
    pepstartExon <- width(prex)[start_ex] + start(peprngCdna) - j[start_ex]
    pependExon <- width(prex)[end_ex] + end(peprngCdna) - j[end_ex]
    ## Transform these coordinates into genomic coords:
    if (strand_num < 0) {
        pependGnm <- end(grObj)[start_ex] - (pepstartExon - 1)
        pepstartGnm <- end(grObj)[end_ex] - (pependExon - 1)
    } else {
        pepstartGnm <- start(grObj)[start_ex] + (pepstartExon - 1)
        pependGnm <- start(grObj)[end_ex] + (pependExon - 1)
    }
    peptides_on_genome <- IRanges(start = pepstartGnm, end = pependGnm)
    chr <- as.character(unique(seqnames((grObj))))
    if (length(chr) > 1)
        stop("Provided exons are on different chromosomes!")

    ## Change the way the mcols is generated. Basically take all gene and
    ## transcript related mcols from the grObj and add additional columns.
    ## Eventually we might opt to rename 'accession' into 'protein_id'.
    want_cols <- c("transcript", "gene", "symbol", "tx_id", "tx_name",
                   "tx_biotype", "gene_name", "gene_id", "gene_biotype")
    got_cols <- colnames(mcols(grObj)) %in% want_cols
    ## Build the mcols:
    mcol <- DataFrame(pepseq = as.character(pfeatures(pObj, pcol = pcol)[[1]]),
                      accession = seqnames(pObj)[1],
                      exonJunctions = junc)
    if (any(got_cols)) {
        to_add <- unique(mcols(grObj)[, got_cols, drop = FALSE])
        if (nrow(to_add) > 1)
            stop("Gene and transcript related mcols of GRanges 'grObj' are not",
                 " unique!")
        mcol <- cbind(to_add, mcol)
    }
    x <- GRanges(seqnames = rep(chr, length(peptides_on_genome)),
                 ranges = peptides_on_genome,
                 strand = strand(grObj)@values,
                 mcol)
    names(x) <- names(peprngProt)
    ## lifting seqinfo over
    seqinfo(x) <- seqinfo(grObj)[chr]

    if (any(mcols(x)$exonJunctions)) {
        rep_num <- (end_ex - start_ex) + 1
        new_x <- pintersect(rep(x, rep_num),
                            grObj[unlist(mapply(start_ex, end_ex, FUN = seq))])
        mcols(new_x)$group <- rep(1:length(x), rep_num)
        ## Get rid of the "hit" col:
        mcols(new_x) <- mcols(new_x)[, colnames(mcols(new_x)) != "hit"]
        return(sort(new_x))
    } else
        mcols(x)$group <- 1:length(x)
    if (validObject(x))
        return(x)
}


## setMethod("mapToGenome", c("Proteins", "GenomicRanges"),
##           function(x, genome, ...) .mapToGenome(x, genome, ...))

setMethod("pmapToGenome", c("Proteins", "GRangesList"),
          function(x, genome, pcol, drop.empty.ranges = TRUE, ...) {
              if (length(x) != length(genome))
                  stop("'x' and 'genome' must have the same length")

              l <- vector("list", length = length(x))
              names(l) <- seqnames(x)
              for (i in seq_len(length(x)))
                  l[[i]] <- tryCatchMapToGenome(x[i], genome[[i]], pcol, ...)
              ans <- GRangesList(l)
              if (drop.empty.ranges)
                  ans <- ans[elementNROWS(ans) > 0]
              if (validObject(ans))
                  return(ans)
          })


## Testing a parallel implementation of the mapping (see issue #)
## benchmarking results are however not that promising (see
## benchmark_pmapToGenome function in test_mapToGenome-ensembldb.R).
setGeneric("pmapToGenome2",
           function(x, genome, ...) standardGeneric("pmapToGenome2"))
setMethod("pmapToGenome2", c("Proteins", "GRangesList"),
          function(x, genome, pcol, drop.empty.ranges = TRUE, ...) {
              if (length(x) != length(genome))
                  stop("'x' and 'genome' must have the same length")
              l <- bpmapply(split(x, 1:length(x)), genome,
                            FUN = tryCatchMapToGenome, ...)
              ## Or better use seqnames(x)?
              names(l) <- seqnames(x)
              ans <- GRangesList(l)
              if (drop.empty.ranges)
                  ans <- ans[elementNROWS(ans) > 0]
              if (validObject(ans))
                  return(ans)
          })


setMethod("mapToGenome", c("Proteins", "GRangesList"),
          function(x, genome, pcol, drop.empty.ranges = TRUE, ...) {
              if (length(x) == 1 & length(genome) == 1) {
                  ans <- tryCatchMapToGenome(x, genome[[1]], ...)
                  if (drop.empty.ranges & length(ans) == 1)
                      ans <- GRangesList()
                  else {
                      ans <- GRangesList(ans)
                      names(ans) <- seqnames(x)
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
                  ans <- pmapToGenome(x = x, genome = genome, pcol = pcol,
                                      drop.empty.ranges = drop.emtpy.ranges, ...)
                  ## get original seqnames order
                  ans <- ans[order(match(names(ans), nmsx))]
              }

              if (validObject(ans))
                  return(ans)

          })


## The method first tries to fetch GRanges representing the CDS encoding the
## proteins from the database and subsequently performs the mapping.
##'
##' @title Map proteins to genomic coordinates using ensembldb
##'
##' @description This method enables the mapping of peptide features within a
##' \code{\linkS4class{Proteins}} object to the genome using annotations
##' provided in an \code{\link[ensembldb]{EnsDb}} object.
##'
##' @details The method tries first to fetch the exons representing the region
##' encoding each of the proteins in \code{x} from the
##' \code{\link[ensembldb]{EnsDb}} provided with parameter \code{genome} and
##' performs then the mapping. To this end, the \code{\linkS4class{Proteins}}
##' object has to provide the necessary IDs that can be used to query the
##' database. Ideally, these should be Ensembl protein IDs.
##'
##' While the mapping between Ensembl protein and transcript IDs is 1:1, a single
##' Uniprot ID can be annotated to several proteins and hence transcripts.
##' If Uniprot IDs are provided, the method thus identifies for 1:n mappings
##' between Uniprot ID and transcript the best matching transcript by comparing
##' the length of the transcripts' CDS with the length of the proteins amino
##' acid sequence.
##'
##' @param x A \code{\linkS4class{Proteins}} object providing the proteins,
##' peptide features to map and IDs to identify the proteins and that can be
##' used to query the database.
##'
##' @param genome A \code{\link[ensembldb]{EnsDb}} object providing the required
##' annotations to perform the mapping.
##'
##' @param id A character vector of length one indicating which metadata columns
##' in \code{x} provide the IDs for the mapping to the transcripts encoding the
##' proteins. Can be the name of any column in \code{acols(x)} or \code{"name"}
##' in which case the \code{seqnames(x)} will be used.
##'
##' @param idType A character vector of length one specifying the type of the
##' provided IDs. Supported are \code{"protein_id"} (Ensembl protein ID),
##' \code{"tx_id"} (Ensembl transcript ID) or \code{"uniprot_id"} (Uniprot ID).
##'
##' @param drop.empty.ranges Wheter empty \code{\link[GenomicFeatures]{GRanges}}
##' should be dropped.
##'
##' @return A \code{\link[GenomicFeatures]{GRangesList}} (names representing the
##' Ensembl transcript IDs) with \code{\link[GenomicFeatures]{GRanges}}
##' representing the genomic coordinates for each peptide feature in \code{x}.
##' The ordering of the results matches the ordering of \code{x}, unless
##' \code{drop.empty.ranges = TRUE} in which case no \code{GRanges} will be
##' returned for proteins that either have no peptide features, or for which
##' the mapping failed. The \code{group} metadata columns of the \code{GRanges}
##' links ranges for peptide features spanning two or more exons.
##'
##' @author Johannes Rainer
##' @examples
##' ## Load the test data
##' data(p)
##'
##' ## Load the EnsDb database for the mapping. In the example human Ensembl
##' ## version 86.
##' library(EnsDb.Hsapiens.v86)
##' edb <- EnsDb.Hsapiens.v86
##'
##' ## The Proteins object contains 9 proteins with each having several peptide
##' ## sequences. The seqnames of the object contain the protein's Uniprot IDs.
##' seqnames(p)
##'
##' ## Below we map these peptide sequences to the genome.
##' res <- mapToGenome(p, edb, idType = "uniprot_id")
##'
##' res
##'
##' ## With the exception of one protein (which Uniprot ID could not be found in
##' the database) for all proteins the mapping could be performed.
##' @noRd
setMethod("mapToGenome", c("Proteins", "EnsDb"),
          function(x, genome, pcol, id = "name", idType = "protein_id",
                   drop.empty.ranges = TRUE, ...) {
              if (!hasProteinData(genome))
                  stop("The submitted database does not provide protein",
                       " annotations.")
              ## Check other input arguments.
              id <- match.arg(id, c("name", avarLabels(x)))
              idType <- match.arg(idType, c("protein_id", "uniprot_id", "tx_id"))

              ## (1) Fetch coding regions for the proteins.
              if (id == "name") {
                  got_id <- seqnames(x)
              } else {
                  got_id <- acols(x)[, id]
              }
              message("Fetch coding region for proteins ... ", appendLF = FALSE)
              cdss <- cdsBy(genome, by = "tx",
                            filter = .featureToFilter(x = idType, got_id),
                            columns = unique(c(idType, "tx_id", "protein_id")))
              message("OK")
              ## Calculate peptide widths - might need them later.
              prot_widths <- width(aa(x))
              names(prot_widths) <- got_id
              tx_widths <- sum(width(cdss))

              if (length(cdss) == 0)
                  stop("Could not find any coding region for the specified IDs.")

              ## (2) Check if we've got results for all input IDs and re-order
              ##     the results
              cdss <- unlist(cdss, use.names = FALSE)
              db_id <- unique(mcols(cdss)[, idType])
              if (!all(got_id %in% db_id))
                  warning("From the ", length(got_id), " provided IDs, ",
                          sum(!(got_id %in% db_id)), " could not be found in ",
                          "database! Peptide features for those can thus not ",
                          "be mapped.")
              common_id <- intersect(got_id, db_id)
              ## Subset the input object and re-order the results.
              keep_what <- which(got_id %in% common_id)
              x <- x[keep_what]
              got_id <- got_id[keep_what]

              cdss <- split(cdss, f = mcols(cdss)[, idType])
              ## How to handle Uniprot: uniprot could be encoded by more than
              ## one tx and each tx can be assigned to more than one uniprot!
              ## Split by uniprot and lapply to return all entries for the
              ## transcript with the best matching coding sequence length.
              if (idType == "uniprot_id") {
                  ## Select for each Uniprot the transcript with the best
                  ## matching CDS length.
                  cdss <- endoapply(cdss, function(z) {
                      diffs <- abs(prot_widths[mcols(z)[1, idType]] -
                                   tx_widths[unique(z$tx_id)]/3)
                      return(z[z$tx_id ==
                               unique(z$tx_id)[which.min(diffs)[1]]])
                  })
              }
              ## Re-order and ensure that the lengths of cdss and x match.
              cdss <- cdss[got_id]

              ## Use parallel processing instead?
              ans <- pmapToGenome(x = x, genome = cdss, pcol = pcol,
                                  drop.empty.ranges = drop.empty.ranges, ...)
              if (validObject(ans))
                  return(ans)
          })



.featureToFilter <- function(x, ...) {
    if (x == "tx_id")
        return(TxidFilter(...))
    if (x == "protein_id")
        return(ProteinidFilter(...))
    if (x == "uniprot_id")
        return(UniprotidFilter(...))
    stop("No filter object for feature ", x)
}
