##' .readAAStringSet
##' Same as Biostrings::readAAStringSet but adds some more metadata to mcols
##' @param filenames fasta filenames
##' @param ... further Arguments passed to Biostrings::readAAStringSet
##' @return AAStringSet
##' @noRd
.readAAStringSet <- function(filenames, ...) {
    isExistingFile <- file.exists(filenames)

    if (any(!isExistingFile)) {
      stop("The file(s) ", paste0(sQuote(filenames[!isExistingFile]),
                                  collapse = ","), " do(es) not exist!")
    }

    ## readAAStringSet can handle multiple input files but we want to know which
    ## entry belongs to which file (and this information is not stored).
    aa <- lapply(filenames, readAAStringSet, ...)
    filenames <- Rle(factor(filenames), lengths = lengths(aa))

    aa <- do.call(c, aa)
    aa <- .addFastaInformation2mcol(aa, fastacomments = names(aa),
                                    filenames = filenames)

    ## we reorder the new Proteins object by the seqnames (AccessionNumber)
    names(aa) <- mcols(aa)$AccessionNumber
    aa <- aa[order(names(aa))]
    aa
}
