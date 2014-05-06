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
  aa <- .addFastaInformation2Mcol(aa, fastacomments = names(aa),
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
                                               name = ametadata(object)$ID[i])
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

#' @param x AAStringSet should be the Proteins@aa slot (subject)
#' @param y MSnExp@featureData data.frame (see
#' fData)
#' @return ProteinCoverageSummary instance
#' @noRd
.proteinCoverageMSnExp <- function(x, y, ..., verbose = TRUE) {
  y <- y[!is.na(y$pepseq), ]

  ## if we have double matches MSnbase combines the entries accession and
  ## description using its utils.vec2ssv function, e.g. "accession1;accession2"
  ## TODO: where is the right place for vec2ssv/ssv2vec ?
  accession <- MSnbase:::utils.ssv2list(y$accession)
  description <- MSnbase:::utils.ssv2vec(y$description)

  m <- vapply(accession, length, integer(1L))
  ## duplicate peptide sequences and databaseFiles
  m <- rep(1:nrow(y), m)
  pepseq <- y$pepseq[m]
  fastafile <- y$databaseFile[m]

  aa <- AAStringSet(pepseq)
  fastacomments <- paste(accession, description)
  aa <- .addMcolAAStringSet(aa, fastacomments = fastacomments,
                            filenames = fastafile)
  .proteinCoverage(x, aa, ..., verbose = verbose)
}

#' @param x AAStringSet should be the Proteins@aa slot (subject)
#' @param y mzID data.frame (see flatten)
#' fData)
#' @return ProteinCoverageSummary instance
#' @noRd
.proteinCoverageMzId <- function(x, y, ..., verbose = TRUE) {
  aa <- AAStringSet(pepseq)
  fastacomments <- paste(y$accession, y$description)
  aa <- .addMcolAAStringSet(aa, fastacomments = fastacomments,
                            filenames = y$databaseFile)
  .proteinCoverage(x, aa, ..., verbose = verbose)
}

#' @param x AAStringSet should be the Proteins@aa slot (subject)
#' @param y AAStringSet or AAStringSetList (pattern)
#' @return ProteinCoverageSummary instance
#' @noRd
.proteinCoverage <- function(x, y, ..., verbose = TRUE) {
  .ProteinCoverageSummary(pattern = y, subject = x, ..., verbose = verbose)
}

