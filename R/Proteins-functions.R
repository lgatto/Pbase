## constructors
.ProteinsFromFasta <- function(filenames, ...) {
  isExistingFile <- file.exists(filenames)

  if (any(!isExistingFile)) {
    stop("The file(s) ", paste0(sQuote(filenames[!isExistingFile]),
                                collapse = ","),
         " do(es) not exist!")
  }

  seqs <- lapply(filenames, readAAStringSet, ...)
  aa <- do.call(c, seqs)

  filenames <- Rle(factor(filenames), lengths = sapply(seqs, length))
  seqnames <- names(seqs)

  nz <- nzchar(seqnames)

  if (!any(nz)) {
    seqnames[!nz] <- as.character(seq(sum(nz)))
  }

  ametadata <- DataFrame(filenames = filenames)
  aa@elementMetadata <- ametadata

  metadata <- list(created = date())

  new("Proteins", aa = aa, metadata = metadata)
}

