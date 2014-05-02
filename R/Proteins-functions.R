## constructors
.ProteinsFromFasta <- function(filenames, ...) {
  isExistingFile <- file.exists(filenames)

  if (any(!isExistingFile)) {
    stop("The file(s) ", paste0(sQuote(filenames[!isExistingFile]),
                                collapse = ","),
         " do(es) not exist!")
  }

  seqs <- lapply(filenames, readAAStringSet, ...)
  seq <- do.call(c, seqs)

  filenames <- Rle(factor(filenames), lengths = sapply(seqs, length))

  mcols <- DataFrame(filenames = filenames)

  metadata <- list(created = date())

  new("Proteins", metadata = metadata, mcols = mcols, seq = seq)
}

