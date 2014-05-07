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

setMethod("[", "Proteins",
          function(x, i, j, ..., drop) {
            p <- callNextMethod()

            if (length(x@pranges)) {
              p@pranges <- x@pranges[i]
            }

            return(p)
          })

## accessor
setMethod("pfeatures",
          "Proteins",
          function(x) extractAt(aa(x), unname(pranges(x))))

setMethod("pranges",
          "Proteins",
          function(x) {
            if (length(x@pranges)) {
              return(x@pranges)
            } else {
              stop("No peptide features found. Do you want to cleave first?")
            }
          })

setMethod("pmetadata",
          "Proteins",
          function(x) lapply(x@pranges, mcols))

setMethod("seqnames",
          signature(x = "Proteins"),
          function(x) mcols(x@aa)$AccessionNumber)

setMethod("show",
          "Proteins",
          function(object) {

  callNextMethod(object)

  cat("Sequences:", tail(capture.output(object@aa), -1), sep = "\n")
})

## replacement

## internal use only; not exported
#setMethod("addpcol",
#          "Proteins",
#          function(x, column, content, force = FALSE) {
#            mcols(x@pranges) <- .addColumn(mcols(x@pranges),
#                                           column = column,
#                                           content = content,
#                                           force = force)
#            x
#          })


## methods
setMethod("cleave",
          "Proteins",
          function(x, enzym = "trypsin", missedCleavages = 0) {
            x@pranges <- cleavageRanges(x = x@aa, enzym = enzym,
                                        missedCleavages = missedCleavages)
            return(x)
          })

setMethod("plot",
          signature(x = "Proteins", y = "missing"),
          function(x, y, ...) .plotProteins(x, ...))

setMethod("proteinCoverage",
          signature(x = "Proteins", y = "mzID"),
          function(x, y, ..., verbose = TRUE) {
            .proteinCoverageMzId(x, flatten(y), ..., verbose = verbose)
          })

