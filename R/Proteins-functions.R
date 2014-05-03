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
  seqnames <- names(aa)

  fastaMetaData <- .fastaComments2DataFrame(seqnames)
  fastaMetaData$filename <- filenames

  ametadata <- fastaMetaData
  aa@elementMetadata <- ametadata

  metadata <- list(created = date())

  new("Proteins", aa = aa, metadata = metadata)
}

.plotProteins <- function(object, from = 1L,
                          to = max(elementLengths(object@aa)), ...) {

  nTracks <- 3L
  tracks <- vector(mode="list", length=length(object) * nTracks)

  for (i in seq(along = object@aa)) {
    idx <- (i - 1L) * nTracks
    tracks[[idx + 1L]] <- ProteinAxisTrack(addNC = TRUE)
    tracks[[idx + 2L]] <- ProteinSequenceTrack(sequence = object@aa[[i]],
                                               name = ametadata(object)$ID[i])
    if (length(object@pfeatures)) {
      ## TODO: adding an ATrack results in an error if "[" is set:
      ## Error in callNextMethod(x, i) :
      ##    bad object found as method (class “function”)
      tracks[[idx + 3L]] <- ATrack(start = start(object@pfeatures[[i]]),
                                   end = end(object@pfeatures[[i]]),
                                   name = "cleavage products")
    }
  }
  tracks <- tracks[sapply(tracks, length) > 0L]

  plotTracks(tracks, from = from, to = to)
}

