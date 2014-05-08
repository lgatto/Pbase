setMethod("Proteins",
          signature(file = "character", uniprotIds = "missing"),
          function(file, uniprotIds, ...) {
              p <- .ProteinsFromFasta(filenames = file, ...)
              names(p@aa) <- acols(p)$AccessionNumber
              p
          })

setMethod("Proteins",
          signature(file = "missing", uniprotIds = "character"),
          function(file, uniprotIds, ...) {
              .toBeImplemented()
          })

## TODO:/BUG: commented because Gviz crashes with callNextMethod error
## see also: https://stat.ethz.ch/pipermail/bioc-devel/2014-May/005701.html
setMethod("[", "Proteins",
          function(x, i, j = "missing", ..., drop) {
              if (!missing(j) || length(list(...)) > 0L)
                  stop("invalid subsetting")
              if (missing(i) || (is.logical(i) && all(i)))
                  return(x)
              if (is.logical(i))
                  i <- which(i)
              if (!is.numeric(i) || any(is.na(i)))
                  stop("invalid subsetting")
              if (any(i < 1) || any(i > length(x)))
                  stop("subscript out of bounds")
              p@aa <- p@aa[i]
              if (length(x@pranges))
                  x@pranges <- x@pranges[i]
              return(x)
          })

## accessor
setMethod("pfeatures", "Proteins",
          function(x) extractAt(aa(x), unname(pranges(x))))

setMethod("pranges", "Proteins",
          function(x) {
            if (length(x@pranges)) {
              return(x@pranges)
            } else {
              stop("No peptide features found. Do you want to cleave first?")
            }
          })


setMethod("length", "Proteins",
          function(x) length(x@aa))

setMethod("metadata", "Proteins",
          function(x) x@metadata)

setMethod("pmetadata", "Proteins",
          function(x) lapply(x@pranges, mcols))

setMethod("ametadata", "Proteins",
          function(x) mcols(x@aa))

setMethod("seqnames","Proteins",
          function(x) names(aa(x)))

setMethod("[[", "Proteins",
          function(x, i, j = missing, ..., drop = TRUE) return(x@aa[[i]]))

setMethod("aa", "Proteins", function(x) x@aa)

## methods
setMethod("addIdentificationData",
          "Proteins",
          function(object, filename) {
            .addIdentificationDataProteins(object, filename)
          })

setMethod("cleave", "Proteins",
          function(x, enzym = "trypsin", missedCleavages = 0) {
            x@pranges <- cleavageRanges(x = x@aa, enzym = enzym,
                                        missedCleavages = missedCleavages)
            x@pranges <- IRangesList(lapply(x@pranges, function(r) {
              mc <- missedCleavages[cumsum(start(r) == 1L)]
              mcols(r) <- DataFrame(MissedCleavages = mc)
              r
            }))
            return(x)
          })

setMethod("plot",
          signature(x = "Proteins", y = "missing"),
          function(x, y, ...) .plotProteins(x, ...))

setMethod("proteinCoverage",
          signature(x = "Proteins", y = "missing"),
          function(x, y, ...) {
            .proteinCoverageProteins(x, ...)
          })

setMethod("show", "Proteins",
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



## internal use only; not exported
setMethod("addacol", "Proteins",
          function(x, column, content, force = FALSE) {
            mcols(x@aa) <- .addColumn(mcols(x@aa),
                                      column = column,
                                      content = content,
                                      force = force)
            x
          })

setMethod("addpcol", "Proteins",
         function(x, column, content, force = FALSE) {
           mcols(x@pranges) <- .addColumn(mcols(x@pranges),
                                          column = column,
                                          content = content,
                                          force = force)
           x
         })

setMethod("aaranges",
          "Proteins",
          function(x, unshift = FALSE) {
            .aarangesProteins(x, unshift = unshift)
          })

