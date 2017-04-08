## Proteins:
## Version 0.1.0 -> 0.2.0
## Change: Remove pranges slot, which are now columns in aa's mcols
##         (see issue #38).


setMethod("updateObject", signature(object = "Proteins"),
          function(object, ..., verbose = TRUE) {
              if (verbose) message("updateObject(object = 'Proteins')")
              object <- asS4(object)
              if (isVersioned(object) && isCurrent(object)["Proteins"])
                callNextMethod()
              else {
                  to <- new(class(object))
                  if (object@.__classVersion__["Proteins"] == "0.1.0") {
                      to@metadata <- object@metadata
                      to@aa <- object@aa
                      mcols(to@aa)$pranges <- pranges(object)
                  }
              }
              if (validObject(to))
                  to
          })
