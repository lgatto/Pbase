##' @param x Proteins object
##' @param unshift if TRUE the IRanges will shift back to start with 1L
##' @return named \code{IRangesList}, \code{length == length(x@@aa)},
##' each element is an \code{IRanges} object starting at 1 and ending
##' at \code{length(x@@aa[i])}.
##' @noRd
.aarangesProteins <- function(x, unshift = FALSE) {
  r <- unname(as(aa(x)@ranges, "IRanges"))
  irl <- split(r, f = seq_along(r))
  names(irl) <- seqnames(x)
  irl
}

##' @param x
##' @param filenames mzIdentML filenames
##' @return a modified Proteins object
##' @noRd
.addIdentificationDataProteins <- function(x, filenames, rmEmptyRanges, par) {
    if (par@IdReader == "mzID") {
        y <- mzID(filenames)
        y <- flatten(y)

        an <- y$accession
        ir <- IRanges(start = y$start, end = y$end)
        names(ir) <- unlist(lapply(strsplit(an, "\\|"), "[", 2))
        fasta <- .fastaComments2DataFrame(paste(y$accession, y$description))
        meta <- as(y[, !colnames(y) %in% c("accession", "description")],
                   "DataFrame")
        mcols(ir) <- cbind(fasta, meta)
        mcols(ir)$filenames <-
                    mcols(ir)$spectrumFile <- Rle(factor(mcols(ir)$spectrumFile))
        mcols(ir)$databaseFile <- Rle(factor(mcols(ir)$databaseFile))
    } else { ## mzR
        .ir <- function(f) {
            if (v) message("  ", k, ". ", f)
            k <<- k + 1
            tmp <- openIDfile(f)
            on.exit(rm(tmp))
            ir <- IRanges()
            if (length(tmp) > 0) {
                y <- psms(tmp)
                an <- as.character(y$DatabaseAccess)
                ir <- IRanges(start = y$start, end = y$end)
                names(ir) <- unlist(lapply(strsplit(an, "\\|"), "[", 2))

                fasta <- .fastaComments2DataFrame(paste(an, y$DatabaseDescription))
                meta <-
                    as(y[, !colnames(y) %in% c("DatabaseAccession", "DatabaseDescription")],
                       "DataFrame")
                mcols(ir) <- cbind(fasta, meta)
            }
            ir
        }
        v <- par@verbose
        k <- 1
        if (v) message("Reading ", length(filenames), " identification files:")
        irl <- lapply(filenames, .ir)
        if (v) message("done.")

        ir <- Reduce(c, irl)
        ir@elementMetadata$filenames <- Rle(factor(filenames),
                                            lengths = lengths(irl))
    }

    nms <- sort(unique(names(ir)))
    mcols(x@aa)$Peptides <- split(ir, names(ir))
    x@aa@elementMetadata$npeps <- elementNROWS(pranges(x))

    if (rmEmptyRanges) {
        x <- rmEmptyRanges(x)
    }
    x
}

##' @param x
##' @param filenames fasta files
##' @return a modified Proteins object
##' @noRd
.addPeptideFragmentsProteins <- function(x, filenames, rmEmptyRanges, par) {
    if (!isEmpty(x@pranges)) {
        stop("The ", sQuote("pranges"), " slot is not empty! ",
            "No ranges and metadata could be added.")
    }
    fragments <- .readAAStringSet(filenames)

    ir <- unlist(.peptidePosition(fragments, x@aa))
    mcols(ir) <- cbind(mcols(fragments)[mcols(ir)$PeptideIndex, ],
                       mcols(ir))
    x@pranges[unique(mcols(ir)$ProteinIndex)] <-
        split(ir, mcols(ir)$ProteinIndex)
    x@aa@elementMetadata$npeps <- elementNROWS(pranges(x))

    if (rmEmptyRanges) {
        x <- rmEmptyRanges(x)
    }
    x
}

##' @param x Proteins object
##' @param mass numeric, length == 2, mass range
##' @param length numeric, length == 2, length range
##' @return modified Proteins object (pcols(x) gains a new "Filtered" column)
##' @noRd
.pfilterProteins <- function(x, mass = NULL, len = NULL) {
    if (isEmpty(x@pranges)) {
        stop("The ", sQuote("pranges"), " slot is empty!")
    }

    filtered <- !.isValidPeptide(pfeatures(x), mass = mass, len = len)
    addpcol(x, "Filtered", unlist(filtered), force = TRUE)
}

.plotProteins <- function(object, from = 1L,
                          to = max(elementNROWS(object@aa)), ...) {

    nTracks <- 3L ## ProteinAxisTrack + ProteinSequenceTrack + 1 prange
    tracks <- vector(mode="list", length=length(object) * nTracks)
    snms <- seqnames(object)

    isRng <- FALSE
    if (ncol(pranges(object)) > 0) {
        isRng <- TRUE
        prngs <- pranges(object)[[1]] ## take first pranges FIXME: use pcol
    }
    
    for (i in seq(along = object@aa)) {
        idx <- (i - 1L) * nTracks
        tracks[[idx + 1L]] <- ProteinAxisTrack(addNC = TRUE,
                                               name = paste0("axis-", snms[i]))

        tracks[[idx + 2L]] <- ProteinSequenceTrack(sequence = object@aa[[i]],
                                                   name = snms[i])

        if (isRng && length(prngs[[i]])) {
            ## TODO: adding an ATrack results in an error if "[" is set:
            ## Error in callNextMethod(x, i) :
            ##    bad object found as method (class “function”)
            tracks[[idx + 3L]] <- ATrack(start = start(prngs[[i]]),
                                         end = end(prngs[[i]]),
                                         name = "peptides",
                                         ...)
        }
    }
    ## ProteinAxisTrack returns length == 0L; that's why we are using the
    ## `is(track[[i]], "NULL")` function here, to exclude empty elements
    tracks <- tracks[!sapply(tracks, is.null)]

    plotTracks(tracks, from = from, to = to)
}

proteotypic <- function(x) {
    stopifnot(inherits(x, "Proteins"))
    if (isEmpty(x@pranges)) {
        stop("The ", sQuote("pranges"), " slot is empty!")
    }

    pp <- as.character(unlist(pfeatures(x)))
    proteotypic <- Rle(pp %in% .singular(pp))
    addpcol(x, "Proteotypic", proteotypic, force = TRUE)
}

rmEmptyRanges <- function(x) {
    x[as.logical(elementNROWS(pranges(x)))]
}

isCleaved <- function(x, missedCleavages = 0) {
    !isEmpty(pranges(x)) &&
        all(missedCleavages %in% unlist(runValue(pmetadata(x)[, "MissedCleavages"])))
}

### Caution: This is based purley on IRanges. No sequence based checks are
### involved! You have to make sure that you compare compareable sequences.
proteinCoverage <- function(x) {
    stopifnot(is(x, "Proteins"))
    prtl <- width(aa(x))
    pepl <- sapply(width(reduce(pranges(x))), sum)
    addacol(x, "Coverage", pepl/prtl)
}


##' @param object An object of class Proteins
##' @param value A new pranges of class CompressedIRangesList
##' @return Proteins object with updated pranges
##' @noRd
replacePranges <- function(object, value) {
    if (length(pranges(object)) != length(value))
        stop("Length of replacement pranges differs from current ones.")
    if (!identical(names(object@pranges), names(value)))
        stop("Names of replacement pranges differ from current ones.")
    object@pranges <- value
    if (validObject(object))
        return(object)
}

##' @param object An object of class Proteins
##' @param value A new acols of class DataFrame
##' @return Proteins object with updated pranges
##' @noRd
replaceAcols <- function(object, value) {
    if (nrow(acols(object)) != nrow(value))
        stop("Number of rows of replacement acols differ from current ones.")
    if (!is.null(rownames(acols(object))) &&
        !identical(rownames(acols(object)), rownames(values)))
        stop("Row names of replacement acols differ from current ones.")
    mcols(object@aa) <- value
    if (validObject(object))
        return(object)
}

