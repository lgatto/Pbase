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
setGeneric("pranges<-", function(object, value) standardGeneric("pranges<-"))
if (is.null(getGeneric("seqnames")))
  setGeneric("seqnames", function(x, ...) standardGeneric("seqnames"))

setGeneric("pvarLabels", function(object, ...) standardGeneric("pvarLabels"))
setGeneric("avarLabels", function(object, ...) standardGeneric("avarLabels"))

## replacement

## methods
if (is.null(getGeneric("addPeptideFragments")))
  setGeneric("addPeptideFragments",
             function(object, filenames, ...) standardGeneric("addPeptideFragments"))
if (is.null(getGeneric("cleave")))
  setGeneric("cleave", function(x, ...) standardGeneric("cleave"))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
if (is.null(getGeneric("pfilter")))
  setGeneric("pfilter", function(x, y, ...) standardGeneric("pfilter"))

## Ranges
setGeneric("proteinCoding",
           function(object, ...) standardGeneric("proteinCoding"))
