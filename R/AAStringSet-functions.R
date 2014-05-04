.addMcolAAStringSet <- function(aa, fastacomments, filenames) {

  if (missing(fastacomments)) {
    fastacomments <- names(aa)
  }

  if (!is(filenames, "Rle")) {
    filenames <- Rle(factor(filenames))
  }

  fastaMetaData <- .fastaComments2DataFrame(fastacomments)
  fastaMetaData$filename <- filenames

  mcols(aa) <- fastaMetaData

  return(aa)
}

.coverageAAStringSet <- function(pattern, subject, verbose = TRUE) {
  pid <- mcols(pattern)$AccessionNumber
  sid <- mcols(subject)$AccessionNumber

  if (is.null(pid)) {
    warning("Could not found any AccessionNumber in ", sQuote("pattern"), "!")
    return(RleList(lapply(elementLengths(subject), function(n) {
      Rle(logical(1L), lengths = n)
    })))
  }

  idx <- match(pid, sid)

  patternList <- unname(split(pattern,
                              factor(idx, levels = seq_along(subject))))

  l <- vector(mode = "list", length = length(subject))

  if (verbose) {
    pb <- txtProgressBar(0L, length(subject), style = 3L)
  }

  for (i in seq(along = subject)) {
    l[[i]] <- .coverageAAString(patternList[[i]], subject[[i]])
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbose) {
    close(pb)
  }
  RleList(l)
}

.coverageOutputTextAAStringSet <- function(object, rlel) {
  il <- which(!as(rlel, "LogicalList"))
  object[il] <- AAStringSet(tolower(object[il]))
  object
}
