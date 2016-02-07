setMethod("Proteins", c("missing", "missing"),
          function(file, uniprotIds, ...){
              aa <- new("AAStringSet")
              aa@elementMetadata <- DataFrame()
              new("Proteins",
                  metadata = list(created = date()),
                  aa = aa,
                  pranges = IRangesList())
          })

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

setMethod("[", "Proteins",
          function(x, i, j = "missing", ..., drop) {
              if (!missing(j) || length(list(...)) > 0L)
                  stop("invalid subsetting")
              if (missing(i) || (is.logical(i) && all(i)))
                  return(x)
              if (is.logical(i))
                  i <- which(i)
              if (is.character(i))
                  i <- match(i, seqnames(x))
              if (!is.numeric(i) || any(is.na(i)))
                  stop("invalid subsetting")
              if (any(i < 1) || any(i > length(x)))
                  stop("subscript out of bounds")

              x@aa <- x@aa[i]
              x@pranges <- x@pranges[i]

              return(x)
          })


## accessors
setMethod("pfeatures", "Proteins",
          function(x) extractAt(aa(x), unname(pranges(x))))

setMethod("pranges", "Proteins",
          function(x) x@pranges)

setMethod("length", "Proteins",
          function(x) length(x@aa))

setMethod("metadata", "Proteins",
          function(x) x@metadata)

## signature is (x, ..., value)
setReplaceMethod("metadata", "Proteins",
                 function(x, name, value) {
                     if (name == "created")
                         stop("Creation date can't be modified.")
                   x@metadata[[name]] <- value
                   return(x)
               })

setMethod("pmetadata", "Proteins",
          function(x) {
              if (!is.null(x@pranges@unlistData@elementMetadata)) {
                  f <- rep.int(seqnames(x), elementLengths(x@pranges))
                  split(x@pranges@unlistData@elementMetadata, f)
              } else {
                  return(NULL)
              }
          })

setMethod("ametadata", "Proteins",
          function(x) mcols(x@aa))

setMethod("seqnames","Proteins",
          function(x) names(aa(x)))


setMethod("avarLabels", "Proteins",
          function(object) names(aa(object)@elementMetadata))

setMethod("pvarLabels", "Proteins",
          function(object) names(pranges(object)@unlistData@elementMetadata@listData))

setMethod("[[", "Proteins",
          function(x, i, j = missing, ..., drop = TRUE) return(x@aa[[i]]))

setMethod("aa", "Proteins", function(x) x@aa)


## Methods
setMethod("addIdentificationData",
          c("Proteins", "character"),
          function(object, id, rmEmptyRanges = TRUE, par = Pparams()) {
              .addIdentificationDataProteins(object, filenames = id,
                                             rmEmptyRanges = rmEmptyRanges,
                                             par = par)
          })

setMethod("addPeptideFragments",
          c("Proteins", "character"),
          function(object, filenames, rmEmptyRanges = TRUE, par = Pparams()) {
              .addPeptideFragmentsProteins(object, filenames = filenames,
                                           rmEmptyRanges = rmEmptyRanges,
                                           par = par)
          })

setMethod("cleave", "Proteins",
          function(x, enzym = "trypsin", missedCleavages = 0, ...) {
              rng <- cleavageRanges(x = x@aa, enzym = enzym,
                                          missedCleavages = missedCleavages,
                                          ...)
              x@pranges <- IRangesList(lapply(rng, function(r) {
                  mc <- missedCleavages[cumsum(start(r) == 1L)]
                  mcols(r) <- DataFrame(MissedCleavages = Rle(mc))
                  r
              }))
              return(x)
          })

setMethod("plot",
          signature(x = "Proteins", y = "missing"),
          function(x, y, ...) .plotProteins(x, ...))

setMethod("pfilter",
          "Proteins",
          function(x, mass = NULL, len = NULL, ...) {
              .pfilterProteins(x, mass = mass, len = len, ...)
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
              if (length(object) > 0) {
                  sn <- seqnames(object)
                  ln <- length(object)
                  cat("Sequences:\n  "); htcat(sn, n = 2)
                  cat("Sequence features:\n  "); htcat(avarLabels(object), n = 2)
                  cat("Peptide features:")
                  if (isEmpty(pranges(object))) cat(" None\n")
                  else {
                      cat("\n  ")
                      htcat(pvarLabels(object), n = 2)
                  }
              }
          })


## internal use only; not exported

setMethod("aaranges",
          "Proteins",
          function(x, ...) {
              .aarangesProteins(x, ...)
          })

