################################################################################
## Proteins
################################################################################

## constructor
if (is.null(getGeneric("Proteins")))
  setGeneric("Proteins", function(file, uniprotIds, ...)
             standardGeneric("Proteins"))

## accessors
if (is.null(getGeneric("aa")))
  setGeneric("aa", function(x, ...) standardGeneric("aa"))
if (is.null(getGeneric("aaranges")))
  setGeneric("aaranges", function(x, ...) standardGeneric("aaranges"))
if (is.null(getGeneric("ametadata")))
  setGeneric("ametadata", function(x, ...) standardGeneric("ametadata"))
if (is.null(getGeneric("acols")))
  setGeneric("acols", function(x, ...) ametadata(x, ...))
if (is.null(getGeneric("pmetadata")))
  setGeneric("pmetadata", function(x, ...) standardGeneric("pmetadata"))
if (is.null(getGeneric("pcols")))
  setGeneric("pcols", function(x, ...) pmetadata(x, ...))
if (is.null(getGeneric("pfeatures")))
  setGeneric("pfeatures", function(x, ...) standardGeneric("pfeatures"))
if (is.null(getGeneric("pranges")))
  setGeneric("pranges", function(x, ...) standardGeneric("pranges"))
if (is.null(getGeneric("seqnames")))
  setGeneric("seqnames", function(x, ...) standardGeneric("seqnames"))

## replacement

## methods
if (is.null(getGeneric("cleave")))
  setGeneric("cleave", function(x, ...) standardGeneric("cleave"))
if (is.null(getGeneric("isCleaved")))
  setGeneric("isCleaved", function(x, ...) standardGeneric("isCleaved"))
if (is.null(getGeneric("plot")))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
if (is.null(getGeneric("pfilter")))
  setGeneric("pfilter", function(x, y, ...) standardGeneric("pfilter"))

setGeneric("pvarLabels", function(object, ...) standardGeneric("pvarLabels"))
setGeneric("avarLabels", function(object, ...) standardGeneric("avarLabels"))

## Ranges
setGeneric("proteinCoding",
           function(object, ...) standardGeneric("proteinCoding"))
