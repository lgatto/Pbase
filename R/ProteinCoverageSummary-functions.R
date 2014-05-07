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

#' @param pattern IRangesList of interest (from mzID, ...)
#' @param subject IRangesList to compare against (mostly from Proteins)
#' @return a LogicaList (TRUE for overlap)
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

  ## find overlap
  subjectForest <- IntervalForest(subject)
  hits <- as.matrix(findOverlaps(pattern, subjectForest, maxgap = 0L,
                                 minoverlap = 1L, type = "within",
                                 select = "all"))
  ## turn Hits into a LogicalList
  subwidth <- unlist(width(subject))

  csw <- (cumsum(c(0L, subwidth)))

  ## shift pattern per corresponding protein length (similar to
  ## AAStringSet@ranges)
  pattern <- .flatIRangesList(pattern[hits[, 1L]],
                              shift = TRUE, shiftBy = csw[hits[, 2L]])

  ## create large logical vector to generate coverage mask in a vectorized way
  l <- Rle(logical(1L), lengths = sum(subwidth))
  l[pattern] <- TRUE

  l <- LogicalList(split(l, rep.int(seq_along(subject), subwidth)))
  names(l) <- names(subject)
  l
}

.calculateProteinCoverageSummary <- function(object) {
  m <- setNames(mean(object@coverage), accessionNumber(object))
  attr(m, "overall") <- mean(m)
  m
}

