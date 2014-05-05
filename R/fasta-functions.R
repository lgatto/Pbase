## UniProt definition: http://www.uniprot.org/help/fasta-headers
## Please note the Accession number is according to
## http://www.uniprot.org/faq/6 the stable ID of the entries.
## In contrast the EntryName could change
.fastaCommentParser <- function(x) {
  ## >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName[ GN=GeneName]
  ## PE=ProteinExistence SV=SequenceVersion
  ## e.g.:
  ## tr|B1XAE9|B1XAE9_ECODH Lipoprotein OS=Escherichia coli (strain K12 / DH10B)
  ## GN=nlpB PE=4 SV=1
  ## (of course, no newlines!)
  ##
  ## short perl explanation:
  ## (?<NAME>...) creates a named group
  ## (?:...) creates a group but doesn't list it in the results
  ## (...)? optional group
  ## \\s matches spaces
  rx <- gregexpr(pattern = paste0(
    ## Database ID, e.g. "tr"
    "^>?(?<DB>[a-z]+)\\|",
    ## AccessionNumber/UniqueIdentifier, e.g "B1XAE9"
    "(?<AN>[A-Z_0-9-]+)\\|",
    ## EntryName, e.g. "B1XAE9_ECODH"
    "(?<EN>[A-Z_0-9]+)\\s+",
    ## optional isoform (could be "Isoform IsoformName of ProteiName")
    "(?:Isoform\\s+(?<IF>.*)\\s+of\\s+)?",
    ## ProteinName, e.g. Lipoprotein
    "(?<PN>.+)\\s+",
    ## OrganismName, e.g. Escherichia coli
    "OS=(?<OS>[^=]+)",
    ## GeneName, e.g. nlpB
    "(?:\\s+GN=(?<GN>[^= ]+))?",
    ## ProteinExistence, e.g. 4
    "(?:\\s+PE=(?<PE>[^= ]+))?",
    ## SequenceVersion, e.g. 1
    "(?:\\s+SV=(?<SV>[^= ]+))?$"), text = x, perl = TRUE)

  x <- do.call(rbind, .perlRxSubstring(x, rx))
  x[!nchar(x)] <- NA_character_
  x
}

## Details:http://www.uniprot.org/manual/accession_numbers
## This regexpr already includes the upcoming extension to 10 characters
## (announced for 07/2014: http://www.uniprot.org/changes/accession_format).
.isUniProtAccessionNumber <- function(x) {
  grepl(paste0("^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]",
               "([A-Z][A-Z0-9]{2}[0-9]){1,2}$"), x)
}

.isPbaseAccessionNumber <- function(x) {
  grepl("^Pb[0-9]+$", x)
}

.isValidAccessionNumber <- function(x) {
  isUniProtAccessionNumber(x) | isPbaseAccessionNumber(x)
}

#' these AccessionNumbers replace missing UniProt AccessionNumbers
#' format: Pb[0-9]+
#' @param x numbers
#' @return character, ID
#' @noRd
.createPbaseAccessionNumbers <- function(x) {
  paste0("Pb", x)
}

.perlRxSubstring <- function(x, rx) {
  mapply(function(xx, rr) {
    start <- attr(rr, "capture.start")
    end <- start + attr(rr, "capture.length") - 1L
    substring(xx, start, end)
  }, xx = x, rr = rx, SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

.fastaComments2DataFrame <- function(x) {
  m <- .fastaCommentParser(x)
  comments <- rep(NA_character_, nrow(m))

  ## any missing AccessionNumber?
  if (anyNA(m[, 2L])) {
    isNA <- which(is.na(m[, 2L]))
    m[isNA, 2L] <- .createPbaseAccessionNumbers(seq_along(isNA))
    comments[isNA] <- gsub("^>", "", x[isNA])
  }

  DataFrame(DB = Rle(factor(m[, 1L])),
            AccessionNumber = m[, 2L],
            EntryName = m[, 3L],
            IsoformName = Rle(m[, 4L]),
            ProteinName = m[, 5L],
            OrganismName = Rle(factor(m[, 6L])),
            GeneName = Rle(factor(m[, 7L])),
            ## Levels of evidence:
            ## http://www.uniprot.org/manual/protein_existence
            ProteinExistence = Rle(factor(m[, 8L],
                                      labels = c("Evidence at protein level",
                                                 "Evidence at transcript level",
                                                 "Inferred from homology",
                                                 "Predicted",
                                                 "Uncertain"),
                                      levels = 1L:5L)),
            SequenceVersion = Rle(m[, 9L]),
            Comment = Rle(comments))
}

