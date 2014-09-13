.Pparams <- setClass("Pparams",
                     slots = c(
                         DbFormat = "character",
                         IdFormat = "character",
                         IdReader = "character",
                         verbose = "logical"))

Pparams <- function(DbFormat = "UniProt",
                    IdFormat = "mzid",
                    IdReader = c("mzR", "mzID"),
                    verbose = TRUE) {
    .Pparams(DbFormat = DbFormat,
             IdReader = match.arg(IdReader),
             IdFormat = IdFormat,
             verbose = verbose)
}


setMethod("show", "Pparams",
          function(object) {
              snms <- slotNames(object)
              vals <- sapply(snms, function(s) slot(object, s))
              txt <- paste(sprintf("%s: %s", snms, vals), collapse = "; ")
              cat(strwrap(txt, exdent = 2), sep = "\n")
          })

