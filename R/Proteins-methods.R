setMethod("Proteins",
          signature(file = "character", uniprotIds = "missing"),
          function(file, uniprotIds, ...) {
            .ProteinsFromFasta(filenames = file, ...)
          })

setMethod("Proteins",
          signature(file = "missing", uniprotIds = "character"),
          function(file, uniprotIds, ...) {
            .toBeImplemented()
          })

## BUG: commented because Gviz crashes with callNextMethod error
#setMethod("[", "Proteins",
#          function(x, i, j, ..., drop)
#          {
#            if (!missing(j) || length(list(...)) > 0L) {
#              stop("invalid subsetting")
#            }
#            if (missing(i) || (is.logical(i) && all(i))) {
#              return(x)
#            }
#            if (is.logical(i)) {
#              i <- which(i)
#            }
#            if (!is.numeric(i) || any(is.na(i))) {
#              stop("invalid subsetting")
#            }
#            if (any(i < 1) || any(i > length(x))) {
#              stop("subscript out of bounds")
#            }
#
#            if (length(pfeatures)) {
#              pfeatures <- x@pfeatures[i]
#            } else {
#              pfeatures <- IRangesList
#            }
#
#            new(class(x),
#                aa = x@aa[i],
#                pfeatures = pfeatures,
#                metadata = x@metadata)
#          })

setMethod("[[",
          "Proteins",
          function(x, i, j, ..., drop = TRUE) {
            return(x@aa[[i]])
          })

setMethod("aa",
          "Proteins",
          function(x) x@aa)

setMethod("cleave",
          "Proteins",
          function(x, enzym = "trypsin", missedCleavages = 0) {
            x@pfeatures <- cleavageRanges(x = x@aa, enzym = enzym,
                                          missedCleavages = missedCleavages)
            return(x)
          })

setMethod("pfeatures",
          "Proteins",
          function(x) x@pfeatures)

setMethod("ametadata",
          "Proteins",
          function(x) mcols(x@aa))

setMethod("length",
          "Proteins",
          function(x) length(x@aa))

setMethod("metadata",
          "Proteins",
          function(x) x@metadata)

setMethod("plot",
          signature(x = "Proteins", y = "missing"),
          function(x, y, ...) .plotProteins(x, ...))

setMethod("pmetadata",
          "Proteins",
          function(x) mcols(x@pfeatures))

setMethod("show",
          "Proteins",
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

  cat("Sequences:", tail(capture.output(object@aa), -1), sep = "\n")
})

