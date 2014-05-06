.toBeImplemented <- function() {
  parentCall <- sys.call(-1L)
  stop(deparse(parentCall), " is not implemented yet!", call.=FALSE)
}

#' @param x IRangesList
#' @param shift if TRUE the IRanges will shift (similar to AAStringSet@ranges)
#' @return an IRanges similar to AAStringSet's @ranges slot (an AAStringSet is
#' just a long "character" vector with ranges)
#' @noRd
.flatIRangesList <- function(x, shift = FALSE) {
  stopifnot(is(x, "IRangesList"))
  if (shift) {
    w <- c(0L, unlist(lapply(end(x), tail, n = 1L)))
    x <- shift(x, head(cumsum(w), -1L))
  }
  unlist(x)
}

#' @param x IRanges
#' @param f defines the grouping (see ?split)
#' @param unshift if TRUE the IRanges will shift back to start with 1L
#' @return a IRangesList
#' @noRd
.splitIRanges <- function(x, f = seq_along(x), unshift = FALSE) {
  stopifnot(is(x, "IRanges"))
  x <- split(unname(x), f)

  if (unshift) {
    w <- c(0L, unlist(lapply(end(x), tail, n = 1L)))
    x <- shift(x, -head(w, -1L))
  }
  x
}

