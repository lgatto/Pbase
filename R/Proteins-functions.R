## constructors
.ProteinsFromFasta <- function(filenames, ...) {
  isExistingFile <- file.exists(filenames)

  if (any(!isExistingFile)) {
    stop("The file(s) ", paste0(sQuote(filenames[!isExistingFile]),
                                collapse = ","), " do(es) not exist!")
  }

  ## readAAStringSet can handle multiple input files but we want to know which
  ## entry belongs to which file (and this information is not stored).
  aa <- lapply(filenames, readAAStringSet, ...)
  filenames <- Rle(factor(filenames), lengths = sapply(aa, length))

  aa <- do.call(c, aa)
  aa <- .addFastaInformation2mcol(aa, fastacomments = names(aa),
                                  filenames = filenames)

  ## we reorder the new Proteins object by the seqnames (AccessionNumber)
  names(aa) <- mcols(aa)$AccessionNumber
  aa <- aa[order(names(aa))]

  ## pranges should have the same order and the same names
  pranges <- IRangesList(replicate(length(aa), IRanges()))
  names(pranges) <- names(aa)

  metadata <- list(created = date())

  p <- new("Proteins", aa = aa, pranges = pranges, metadata = metadata)
}

#' @param x Proteins object
#' @param unshift if TRUE the IRanges will shift back to start with 1L
#' @return named IRangesList, length == length(x@aa), each element is an IRanges
#' object starting at 1 and ending at length(x@aa[i]).
#' @noRd
.aarangesProteins <- function(x, unshift = FALSE) {
  r <- unname(as(aa(x)@ranges, "IRanges"))
  irl <- .splitIRanges(r, unshift = unshift)
  names(irl) <- seqnames(x)
  irl
}

#' @param x 
#' @param filename mzIdentML filename
#' @return a modified Proteins object
#' @noRd
.addIdentificationDataProteins <- function(x, filename, rmEmptyRanges) {
  if (!isEmpty(x@pranges)) {
    stop("The ", sQuote("pranges"), " slot is not empty! ",
         "No ranges and metadata could be added.")
  }

  ## y <- flatten(mzID(filename))
  load("~/tmp/y.rda")

  an <- y$accession
  ir <- IRanges(start = y$start, end = y$end)
  names(ir) <- unlist(lapply(strsplit(an, "\\|"), "[", 2))

  fasta <- Pbase:::.fastaComments2DataFrame(paste(y$accession, y$description))
  meta <- as(y[, !colnames(y) %in% c("accession", "description")],
             "DataFrame")
  mcols(ir) <- cbind(fasta, meta)

  if (rmEmptyRanges) {
      x@pranges <- split(ir, names(ir))
      nms <- names(x@pranges)
      x@aa <- aa(x)[nms]
  } else {
      n <- length(x)
      .irl <- IRangesList(replicate(n, IRanges()))
      names(.irl) <- seqnames(x)
      .irl[names(ir)] <- split(ir, names(ir))     
      x@pranges <- .irl
  }
  x
}

#' @param x Proteins object
#' @param mass numeric, length == 2, mass range
#' @param length numeric, length == 2, length range
#' @return modified Proteins object (pcols(x) gains a new "Filtered" column)
#' @noRd
.pfilterProteins <- function(x, mass = NULL, len = NULL) {
    if (isEmpty(x@pranges)) {
        stop("The ", sQuote("pranges"), " slot is empty!")
    }

    filtered <- !.isValidPeptide(pfeatures(x), mass = mass, len = len)
    addpcol(x, "Filtered", unlist(filtered), force = TRUE)
}

.plotProteins <- function(object, from = 1L,
                          to = max(elementLengths(object@aa)), ...) {

  nTracks <- 3L
  tracks <- vector(mode="list", length=length(object) * nTracks)
  snms <- seqnames(object)
  
  for (i in seq(along = object@aa)) {
    idx <- (i - 1L) * nTracks
    tracks[[idx + 1L]] <- ProteinAxisTrack(addNC = TRUE,
                                           name = paste0("axis-", snms[i]))
  
    tracks[[idx + 2L]] <- ProteinSequenceTrack(sequence = object@aa[[i]],
                                               name = snms[i])

    if (length(object@pranges[[i]])) {
      ## TODO: adding an ATrack results in an error if "[" is set:
      ## Error in callNextMethod(x, i) :
      ##    bad object found as method (class “function”)
      tracks[[idx + 3L]] <- ATrack(start = start(object@pranges[[i]]),
                                   end = end(object@pranges[[i]]),
                                   name = "peptides")
    }
  }
  ## ProteinAxisTrack returns length == 0L; that's why we are using the
  ## `is(track[[i]], "NULL")` function here, to exclude empty elements
  tracks <- tracks[!sapply(tracks, is.null)]

  plotTracks(tracks, from = from, to = to)
}

#' @param x Proteins object
#' @return a modified Proteins object
#' @noRd
.proteincoverageproteins <- function(x, ...) {
  .proteinCoverageProteinsRanges(x, ranges = pranges(x), ...)
}

#' @param x Proteins object
#' @param ranges IRanges object containing the ranges of the peptides
#' fData)
#' @return a modified Proteins object
#' @noRd
.proteinCoverageProteinsRanges <- function(x, ranges, ...) {
    subject <- aaranges(x, unshift = TRUE)
    coverage <- .proteinCoverage(pattern = ranges, subject = subject)

    addacol(x, "Coverage", coverage)
}

#' @param x Proteins object
#' @return a modified Proteins object (pcols(x) gains a "Proteotypic" column)
#' @noRd
.proteotypicProteins <- function(x) {
    if (isEmpty(x@pranges)) {
        stop("The ", sQuote("pranges"), " slot is empty!")
    }

    pp <- as.character(unlist(pfeatures(x)))
    proteotypic <- Rle(pp %in% .singular(pp))
    addpcol(x, "Proteotypic", proteotypic, force = TRUE)
}


## might become a method
rmEmptyRanges <- function(x) {
    lns <- elementLengths(pranges(x)); 
    em <- lns == 0
    x[!em]
}
