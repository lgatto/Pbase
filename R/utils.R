.toBeImplemented <- function() {
  parentCall <- sys.call(-1L)
  stop(deparse(parentCall), " is not implemented yet!", call.=FALSE)
}

##' @param x data.frame
##' @param column name of the column
##' @param content content to insert into the new x[, column]
##' @param force if TRUE it overwrites an existing column otherwise an
##' error is thrown
##' @return x with the new column "column"
##' @noRd
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

##' test whether a given vector is in the range from:to
##' @param x numeric
##' @param range numeric, lengh == 2
##' @return logical of length(x)
##' @noRd
.isInRange <- function(x, range) {
    stopifnot(is.numeric(range) && length(range) == 2L)
    if (range[1L] > range[2L]) {
        range[2L:1L] <- range[1L:2L]
    }
    range[1L] <= x & x <= range[2L]
}

##' @param x IRangesList
##' @param shift if TRUE the IRanges will shift (similar to
##' AAStringSet@@ranges)
##' @param shiftBy if shift TRUE and shiftBy it is used to shift the IRanges,
##' lenght must be equal to 1 or the length of x
##' @return an IRanges similar to AAStringSet's ranges slot (an AAStringSet is
##' just a long "character" vector with ranges)
##' @noRd
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

##' setNames2 is similar to setNames but uses the names of the second
##' argument if they are available
##' @param object object to give names to
##' @param nm names/object with names
##' @return named object
##' @noRd
.setNames2 <- function(object, nm) {
    if (is.null(names(nm))) {
        names(object) <- as.character(nm)
    } else {
        names(object) <- names(nm)
    }
    object
}

##' Splits a sequence into single AA
##' @param x character, AAString, AAStringSet
##' @return list of single characters
##' @noRd
.singleAA <- function(x) {
    if (is(x, "AAStringSetList")) {
        x <- unlist(x)
    }
    ## if we always apply as.character the names of original characters are
    ## removed
    if (!is.character(x)) {
        x <- as.character(x)
    }
    strsplit(x, split = character(), fixed = TRUE)
}

##' In contrast to base::unique this function removes all duplicated
##' values and keeps only singular ones. base::unique keeps one copy
##' of each duplicate.
##' @param x vector
##' @return reduced x
##' @noRd
.singular <- function(x) {
    x[! (duplicated(x) | duplicated(x, fromLast = TRUE))]
}

addacol <- function(x, column, content, force = FALSE) {
    mcols(x@aa) <- .addColumn(mcols(x@aa),
                              column = column,
                              content = content,
                              force = force)
    x
}

addpcol <- function(x, column, content, force = FALSE) {
    x@pranges@unlistData@elementMetadata <-
        .addColumn(x@pranges@unlistData@elementMetadata,
                   column = column,
                   content = content,
                   force = force)
    x
}

htcat <- function(x, n = 3) {
    nx <- length(x)
    if (nx <= n) {
        cat(paste(x, collapse = ", "))
    } else {
        cat(paste(paste0("[", 1:n, "]"), head(x, n)))
        cat(" ... ")
        cat(paste(paste0("[", (nx-n+1):nx, "]"), tail(x, n)))
    }
    cat("\n")
}


##' @title Are all the ranges on the same strand
##' @param gr A \code{GRanges} object.
##' @return A logical if \emph{all} the ranges in the \code{gr} object
##' are on the \code{"-"} (or \code{"+"} for code{isForward}) strand.
##' @author Laurent Gatto
##' @aliases isForward
isReverse <- function(gr)
    isTRUE(all(strand(gr) == "-"))

##' @rdname isReverse
isForward <- function(gr)
    isTRUE(all(strand(gr) == "+"))



### modified, from the IRanges vignette
### internal, for testing
plotRanges <- function(x, xlim = x, main = "",
                       col = "black", sep = 0.5, ...) {
    if (is(x, "GRanges"))
        x <- ranges(x)
    height <- 1
    if (is(xlim, "Ranges"))
        xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
    title(main)
    axis(1)
}
