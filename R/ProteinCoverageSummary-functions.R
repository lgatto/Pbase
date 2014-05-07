## constructors
#' @param x Proteins object
#' @param ranges IRanges object containing the ranges of the peptides
#' @noRd
.ProteinCoverageSummary <- function(x, ranges, verbose = TRUE) {

  metadata <- list(created = date())

  subject <- .aaRanges(x, unshift = TRUE)

  coverage <- .calculateProteinCoverage(pattern = ranges, subject = subject)

  new("ProteinCoverageSummary", aa = x@aa, coverage = coverage,
      metadata = metadata)
}

#.areRangesCorrect <- function(aa, ranges, bb) {
#  a <- rep(aa, elementLengths(ranges))
#  all(substring(rep(aa, elementLengths(ranges)), unlist(start(ranges)), unlist(end(ranges))) == bb)
#}

.calculateProteinCoverageByRegexpr <- function(pattern, subject) {
  if (anyDuplicated(names(subject))) {
    stop("No duplicated names for ", sQuote("subject"), " allowed!")
  }

  lpat <- split(pattern, names(pattern))
  ## Caution: this will maybe fail for large AA
  lpat <- lapply(lpat, function(x)paste0(unlist(x), collapse="|"))
  lsub <- split(subject, names(subject))

  m <- match(names(lsub), names(lpat))

  r <- mapply(function(p, s, j) {
    if (!is.null(p)) {
      rx <- gregexpr(p, s)[[1L]]
      i <- which(rx > 0L)
      l <- attr(rx, "match.length")[i]
      rx <- rx[i]
      if (length(rx)) {
        return(matrix(c(rx, rx+l-1, rep.int(j, length(rx))), ncol=3))
      }
    }
  }, p = lpat[m], s = lsub, j = m, SIMPLIFY = FALSE, USE.NAMES = FALSE)

  if (length(r)) {
    r <- do.call(rbind, r)
    ir <- .splitIRanges(IRanges(start = r[, 1L], end = r[, 2L]), f=r[,3L])
    names(ir) <- names(lpat)[unique(r[, 3L])]
  } else {
    ir <- IRanges()
  }
  ir
}

#' Calculates the coverage of patterns in subjects based on the width of the
#' patterns. This method assumes that all patterns lie within a subject and the
#' the names of patterns and subjects containing the correct accession.
#' Caution: This is based purley on IRanges. No sequence based checks are
#' involved! You have to make sure that you compare compareable sequences.
#' @param pattern IRangesList of interest (from mzID, ...)
#' @param subject IRangesList to compare against (mostly from Proteins)
#' @return double, NA for patterns that are outside of subjects.
#' @noRd
.calculateProteinCoverage <- function(pattern, subject) {
  stopifnot(is(pattern, "IRangesList"))
  stopifnot(is(subject, "IRangesList"))

  if (is.null(names(pattern))) {
    stop("No names for ", sQuote("pattern"), " available!")
  } else if (anyDuplicated(names(pattern))) {
    stop("No duplicated names for ", sQuote("pattern"), " allowed!")
  }
  if (is.null(names(subject))) {
    stop("No names for ", sQuote("subject"), " available!")
  } else if (anyDuplicated(names(subject))) {
    stop("No duplicated names for ", sQuote("subject"), " allowed!")
  }

  ## combine overlapping ranges
  pattern <- reduce(pattern)

  ## reorder patterns by subject names (introduces NULL for non matches)
  orderedPattern <- as.list(pattern)[names(subject)]

  ## replace NULL by empty IRanges()
  orderedPattern[unlist(lapply(orderedPattern, is.null))] <- IRanges()
  orderedPattern <- IRangesList(orderedPattern)

  ## test that all patterns lie within the subject
  patternEnd <- as.integer(lapply(end(orderedPattern), tail, n = 1L))
  subjectEnd <- as.integer(lapply(end(subject), tail, n = 1L))
  isOutside <- which(patternEnd > subjectEnd)

  ## calculate coverage
  coverage <- as.double(sum(width(orderedPattern))/width(subject))
  coverage[isOutside] <- NA
  names(coverage) <- names(subject)

  coverage
}

