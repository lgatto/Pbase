.toBeImplemented <- function() {
  parentCall <- sys.call(-1L)
  stop(deparse(parentCall), " is not implemented yet!", call.=FALSE)
}

#' @param x IRangesList
#' @param shift if TRUE the IRanges will shift (similar to AAStringSet@ranges)
#' @param shiftBy if shift TRUE and shiftBy it is used to shift the IRanges,
#' lenght must be equal to 1 or the length of x
#' @return an IRanges similar to AAStringSet's @ranges slot (an AAStringSet is
#' just a long "character" vector with ranges)
#' @noRd
.flatIRangesList <- function(x, shift = FALSE, shiftBy) {
  stopifnot(is(x, "IRangesList"))
  if (shift) {
    if (missing(shiftBy)) {
      shiftBy <- head(c(0L, unlist(lapply(end(x), tail, n = 1L))), -1L)
    }
    x <- shift(x, shiftBy)
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

