.toBeImplemented <- function() {
  parentCall <- sys.call(-1L)
  stop(deparse(parentCall), " is not implemented yet!", call.=FALSE)
}

#' @param x data.frame
#' @param column name of the column
#' @param content content to insert into the new x[, column]
#' @param force if TRUE it overwrites an existing column otherwise an error is
#' thrown
#' @return x with the new column "column"
#' @noRd
.addColumn <- function(x, column, content, force = FALSE) {
  stopifnot(inherits(x, "data.frame") || inherits(x, "DataFrame"))
  stopifnot(is.character(column) && nchar(column) > 0L)
  stopifnot(length(content) > 0L)

  if (column %in% colnames(x) && !force) {
    stop("The column ", sQuote(column), " already exists. Use ",
         sQuote("force = TRUE"), " to overwrite it.")
  }

  x[, column] <- content
  x
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

#' setNames2 is similar to setNames but uses the names of the applied second
#' argument if they are available
#' @param object object to give names to
#' @param nm names/object with names
#' @return named object
.setNames2 <- function(object, nm) {
    if (is.null(names(nm))) {
        names(object) <- as.character(nm)
    } else {
        names(object) <- names(nm)
    }
    object
}
