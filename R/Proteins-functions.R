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

.addIdentificationDataProteins <- function(object, filename) {
  if (!isEmpty(x@pranges)) {
    stop("The ", sQuote("pranges"), " slot is not empty! ",
         "No ranges and metadata could be added.")
  }

  y <- flatten(mzID(filename))

  an <- y$accession
  ir <- IRanges(start = y$start, end = y$end)
  names(ir) <- unlist(lapply(strsplit(an, "\\|"), "[", 2))

  fasta <- .fastaComments2DataFrame(paste(y$accession, y$description))
  meta <- as(y[, !colnames(y) %in% c("accession", "description")],
             "DataFrame")
  mcols(ir) <- cbind(fasta, meta)

  x@pranges <- split(ir, names(ir))

  x
}

.plotProteins <- function(object, from = 1L,
                          to = max(elementLengths(object@aa)), ...) {

  nTracks <- 3L
  tracks <- vector(mode="list", length=length(object) * nTracks)

  for (i in seq(along = object@aa)) {
    idx <- (i - 1L) * nTracks
    tracks[[idx + 1L]] <- ProteinAxisTrack(addNC = TRUE,
                                           name = paste0("axis-",
                                                         seqnames(object[i])))
    tracks[[idx + 2L]] <- ProteinSequenceTrack(sequence = object@aa[[i]],
                                               name = seqnames(object)[i])
    if (length(object@pranges[[i]])) {
      ## TODO: adding an ATrack results in an error if "[" is set:
      ## Error in callNextMethod(x, i) :
      ##    bad object found as method (class “function”)
      tracks[[idx + 3L]] <- ATrack(start = start(object@pranges[[i]]),
                                   end = end(object@pranges[[i]]),
                                   name = "cleavage products")
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
.proteinCoverageProteins <- function(x, ...) {
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

