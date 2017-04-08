################################################################################
## class definitions
################################################################################

## See issue 38 for a discussion on class update from 0.1 to 0.2.
setClass("Proteins",
         contains = "Versioned",
         slots = list(metadata = "list",
                      aa = "AAStringSet"),
         prototype = prototype(
             new("Versioned",
                 versions = c(P = "0.2"))))


################################################################################
## validity checks
################################################################################

.validProteins <- function(object) {
    ## TODO: extend validity checks and write unit tests
    msg <- validMsg(NULL, isCurrent(object))
    if (is.null(msg)) {
        return(TRUE)
    }
    msg
}

setValidity("Proteins", .validProteins)
