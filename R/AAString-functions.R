.coverageAAString <- function(pattern, subject) {
  pattern <- as.character(pattern)
  subject <- as.character(subject)

  r <- Rle(logical(1L) , lengths = nchar(subject))
  n <- nchar(pattern)

  for (i in seq(along = pattern)) {
    rx <- gregexpr2(pattern = pattern[i], text = subject)[[1L]]
    rx <- as.vector(rx[rx > 0L])
    r[IRanges(rx, width=n[i])] <- TRUE
  }
  r
}
