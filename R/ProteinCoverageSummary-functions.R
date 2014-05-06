## constructors
#' @param x Proteins object
#' @param ranges IRanges object containing the ranges of the peptides
#' @noRd
.ProteinCoverageSummary <- function(x, ranges, verbose = TRUE) {

  metadata <- list(created = date())

  new("ProteinCoverageSummary", aa = x@aa, coverage = coverage,
      metadata = metadata)
}

.calculateProteinCoverageSummary <- function(object) {
  il <- as(object@coverage, "LogicalList")
  m <- setNames(mean(il), accessionNumber(object))
  attr(m, "overall") <- mean(m)
  m
}

