## TODO:/BUG: commented because Gviz crashes with callNextMethod error
## see also: https://stat.ethz.ch/pipermail/bioc-devel/2014-May/005701.html
setMethod("[", "VirtualProteins",
          function(x, i, j, ..., drop) {
            if (!missing(j) || length(list(...)) > 0L) {
              stop("invalid subsetting")
            }
            if (missing(i) || (is.logical(i) && all(i))) {
              return(x)
            }
            if (is.logical(i)) {
              i <- which(i)
            }
            if (!is.numeric(i) || any(is.na(i))) {
              stop("invalid subsetting")
            }
            if (any(i < 1) || any(i > length(x))) {
              stop("subscript out of bounds")
            }

            new(class(x),
                aa = x@aa[i],
                metadata = x@metadata)
          })

setMethod("[[",
          "VirtualProteins",
          function(x, i, j, ..., drop = TRUE) {
            return(x@aa[[i]])
          })

setMethod("aa",
          "VirtualProteins",
          function(x) x@aa)

setMethod("accessionNumber",
          "VirtualProteins",
          function(x) {
            mcols(x@aa)$AccessionNumber
          })

setMethod("ametadata",
          "VirtualProteins",
          function(x) mcols(x@aa))

setMethod("length",
          "VirtualProteins",
          function(x) length(x@aa))

setMethod("metadata",
          "VirtualProteins",
          function(x) x@metadata)

setMethod("show",
          "VirtualProteins",
          function(object) {

  topics <- c("S4 class type",
              "Class version",
              "Created",
              "Number of Proteins")

  topics <- format(topics, justify = "left")

  n <- length(object)

  values <- c(class(object),
              tail(as(classVersion(object), "character"), 1L),
              object@metadata$created,
              n)

  cat(paste0(topics, ": ",  values, collapse = "\n"), sep = "\n")
})


## replacement

## internal use only; not exported
setMethod("addacol",
          "VirtualProteins",
          function(x, column, content, force = FALSE) {
            mcols(x@aa) <- .addColumn(x,
                                      column = column,
                                      content = content,
                                      force = force)
            x
          })

