## constructors
## patterns: AAStringSetList to compare against
## subject: Proteins@aa (AAStringSet)
.ProteinCoverageSummary <- function(pattern, subject, verbose = TRUE) {
  coverage <- .coverageAAStringSet(pattern = pattern, subject = subject,
                                   verbose = verbose)

  metadata <- list(created = date())

  new("ProteinCoverageSummary", aa = subject, coverage = coverage,
      metadata = metadata)
}

.calculateProteinCoverageSummary <- function(object) {
  il <- as(object@coverage, "LogicalList")
  m <- setNames(mean(il), accessionNumber(object))
  attr(m, "overall") <- mean(m)
  m
}

