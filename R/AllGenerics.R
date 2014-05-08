################################################################################
## VirtualProteins
################################################################################

## accessors
if (is.null(getGeneric("aa")))
  setGeneric("aa", function(x, ...) standardGeneric("aa"))
if (is.null(getGeneric("aaranges")))
  setGeneric("aaranges", function(x, ...) standardGeneric("aaranges"))
if (is.null(getGeneric("accessionNumber")))
  setGeneric("accessionNumber", function(x, ...)
             standardGeneric("accessionNumber"))
if (is.null(getGeneric("ametadata")))
  setGeneric("ametadata", function(x, ...) standardGeneric("ametadata"))
if (is.null(getGeneric("acols")))
  setGeneric("acols", function(x, ...) ametadata(x, ...))

## replacement
if (is.null(getGeneric("addacol")))
  setGeneric("addacol", function(x, ...) addacol(x, ...))


################################################################################
## Proteins
################################################################################
if (is.null(getGeneric("Proteins")))
  setGeneric("Proteins", function(file, uniprotIds, ...)
             standardGeneric("Proteins"))
if (is.null(getGeneric("proteinCoverage")))
  setGeneric("proteinCoverage", function(x, y, ...)
             standardGeneric("proteinCoverage"))

## accessors
if (is.null(getGeneric("pmetadata")))
  setGeneric("pmetadata", function(x, ...) standardGeneric("pmetadata"))
if (is.null(getGeneric("pcols")))
  setGeneric("pcols", function(x, ...) pmetadata(x, ...))
if (is.null(getGeneric("pfeatures")))
  setGeneric("pfeatures", function(x, ...) standardGeneric("pfeatures"))
if (is.null(getGeneric("pranges")))
  setGeneric("pranges", function(x, ...) standardGeneric("pranges"))

## replacement
#if (is.null(getGeneric("addpcol")))
#  setGeneric("addpcol", function(x, ...) addpcol(x, ...))
if (is.null(getGeneric("addIdentificationData")))
  setGeneric("addIdentificationData", function(object, ...)
             addIdentificationData(object, ...))

## methods
if (is.null(getGeneric("cleave")))
  setGeneric("cleave", function(x, ...) standardGeneric("cleave"))
if (is.null(getGeneric("plot")))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

