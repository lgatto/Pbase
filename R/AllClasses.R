################################################################################
## class definitions
################################################################################

setClass("Proteins",
         contains = c("Versioned"),
         slots = list(
             metadata = "list",
             aa = "AAStringSet",
             pranges = "CompressedIRangesList"),
         prototype = prototype(
          new("Versioned",
              versions = c(Proteins = "0.1"))))

################################################################################
## validity checks
################################################################################

.validProteins <- function(object) {
    ## TODO: extend validity checks and write unit tests
    msg <- validMsg(NULL, isCurrent(object))

    if (anyDuplicated(seqnames(object))) {
        msg <- validMsg(msg, paste0("Duplicated names in ", sQuote("seqnames"),
                                    " are not allowed!"))
    }

    if (length(object@pranges) != length(object@aa)) {
        msg <- validMsg(msg, paste0("Number of IRanges in @pranges and the ",
                                    "number of sequences in @aa do not match!"))
    }

    if (any(names(object@pranges) != names(object@aa))) {
        msg <- validMsg(msg, paste0("Names of IRanges in @pranges and the ",
                                    "names of sequences in @aa do not match!"))
    }

    if (length(object@aa) != nrow(object@aa@elementMetadata)) {
        msg <- validMsg(msg, paste0("Number of rows in the ametadata and the ",
                                    "number of sequences in @aa do not match!"))
    }

    if (!is.null(object@pranges@unlistData@elementMetadata) &&
        nrow(object@pranges@unlistData@elementMetadata) !=
        sum(elementNROWS(object@pranges))) {
        msg <- validMsg(msg, paste0("Number of rows in the pmetadata and the ",
                                    "number of IRanges in @pranges do not ",
                                    "match!"))
    }

    if (is.null(msg)) {
        return(TRUE)
    }

    msg
}

setValidity("Proteins", .validProteins)
