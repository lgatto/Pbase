################################################################################
## Proteins
################################################################################

## constructor
if (is.null(getGeneric("Proteins")))
  setGeneric("Proteins", function(file, uniprotIds, ...)
             standardGeneric("Proteins"))

if (is.null(getGeneric("seqnames")))
  setGeneric("seqnames", function(x, ...) standardGeneric("seqnames"))

## replacement

## Building fails without this one, but unclear why!?!
setGeneric("acols<-", function(object, value) standardGeneric("acols<-"))


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
