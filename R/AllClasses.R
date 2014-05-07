################################################################################
## class definitions
################################################################################

setClass("VirtualProteins",
         contains = c("Versioned", "VIRTUAL"),
         slots = list(
          metadata = "list",
          aa = "AAStringSet"),
         prototype = prototype(
          new("Versioned",
              versions = c(VirtualProteins = "0.1"))))

setClass("Proteins",
         contains = c("Versioned", "VirtualProteins"),
         slots = list(
          pranges = "CompressedIRangesList"),
         prototype = prototype(
          new("Versioned",
              versions = c(Proteins = "0.1"))))

################################################################################
## validity checks
################################################################################

.validProteins <- function(object) {
  ## TODO:
  msg <- validMsg(NULL, isCurrent(object))

  #n <- length(object@seq)

  #if (nrow(object@mcols) != n) {
  #  msg <- validMsg(msg, paste0("Number of rows in the metadata and ",
  #                              "the number of sequences do not match!"))
  #}

  #npfeatures <- length(object@pfeatures)

  #if (npfeatures && npfeatures != n) {
  #  msg <- validMsg(msg, paste0("Number of IRanges in pfeatures and ",
  #                              "the number of sequences do not match!"))
  #}
  if (is.null(msg)) {
    return(TRUE)
  }

  msg
}

setValidity("Proteins", .validProteins)
