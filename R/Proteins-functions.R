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

  metadata <- list(created = date())

  new("Proteins", aa = aa, metadata = metadata)
}

.plotProteins <- function(object, from = 1L,
                          to = max(elementLengths(object@aa)), ...) {

  nTracks <- 3L
  tracks <- vector(mode="list", length=length(object) * nTracks)

  for (i in seq(along = object@aa)) {
    idx <- (i - 1L) * nTracks
    tracks[[idx + 1L]] <- ProteinAxisTrack(addNC = TRUE)
    tracks[[idx + 2L]] <- ProteinSequenceTrack(sequence = object@aa[[i]],
                                               name = ametadata(object)$AccessionNumber[i])
    if (length(object@pranges)) {
      ## TODO: adding an ATrack results in an error if "[" is set:
      ## Error in callNextMethod(x, i) :
      ##    bad object found as method (class “function”)
      tracks[[idx + 3L]] <- ATrack(start = start(object@pranges[[i]]),
                                   end = end(object@pranges[[i]]),
                                   name = "cleavage products")
    }
  }
  tracks <- tracks[sapply(tracks, length) > 0L]

  plotTracks(tracks, from = from, to = to)
}

#' @param x Proteins object
#' @param y mzID data.frame (see flatten)
#' fData)
#' @return a modified Proteins object
#' @noRd
.proteinCoverageProteinsMzId <- function(x, y, ...) {
  an <- y$accession
  ir <- IRanges(start = y$start, end = y$end)
  names(ir) <- unlist(lapply(strsplit(an, "\\|"), "[", 2))

  fasta <- .fastaComments2DataFrame(paste(y$accession, y$description))
  meta <- as(y[, !colnames(y) %in% c("accession", "description")],
             "DataFrame")
  mcols(ir) <- cbind(fasta, meta)

  irl <- split(ir, names(ir))

  if (length(x@pranges)) {
    warning("The ", sQuote("pranges"), " slot is not empty! ",
            "No ranges and metadata added.")
  } else {
    x@pranges <- irl
  }

  .proteinCoverageProteinsRanges(x, ranges = irl, ...)
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

