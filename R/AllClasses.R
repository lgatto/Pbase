################################################################################
## class definitions
################################################################################
setClass("Proteins",
         contains = "Versioned",
         slots = list(
          metadata = "list",   # global metadata
          mcols = "DataFrame",  # element-wise (AAString) metadata
          pfeatures = "IRanges",
          seq = "AAStringSet"),
         prototype = prototype(
          new("Versioned",
              versions = c(Proteins = "0.1"))))

################################################################################
## validity checks
################################################################################

.validProteins <- function(object) {
  msg <- validMsg(NULL, isCurrent(object))

  if (nrow(object@mcols) != length(object@seq)) {
    msg <- validMsg(msg, paste0("Number of rows in the metadata and ",
                                "the length of the sequences do not match!"))
  }

  if (is.null(msg)) {
    return(TRUE)
  }

  msg
}

setValidity("Proteins", .validProteins)
