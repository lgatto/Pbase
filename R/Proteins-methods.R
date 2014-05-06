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

setMethod("cleave",
          "Proteins",
          function(x, enzym = "trypsin", missedCleavages = 0) {
            x@pranges <- cleavageRanges(x = x@aa, enzym = enzym,
                                        missedCleavages = missedCleavages)
            return(x)
          })

setMethod("pfeatures",
          "Proteins",
          function(x) {
            if (length(x@pranges)) {
              return(extractAt(x@aa, x@pranges))
            } else {
              stop("No peptide features found. Do you want to cleave first?")
            }
          })

setMethod("plot",
          signature(x = "Proteins", y = "missing"),
          function(x, y, ...) .plotProteins(x, ...))

setMethod("pmetadata",
          "Proteins",
          function(x) mcols(x@pranges))

setMethod("proteinCoverage",
          signature(x = "Proteins", y = "mzID"),
          function(x, y, ..., verbose = TRUE) {
            .proteinCoverageMzId(x, flatten(y), ..., verbose = verbose)
          })

setMethod("seqnames",
          signature(x = "Proteins"),
          function(x) mcols(x@aa)$AccessionNumber)

setMethod("show",
          "Proteins",
          function(object) {

  callNextMethod(object)

  cat("Sequences:", tail(capture.output(object@aa), -1), sep = "\n")
})

