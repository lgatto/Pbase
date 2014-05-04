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

.perlRxSubstring <- function(x, rx) {
  mapply(function(xx, rr) {
    start <- attr(rr, "capture.start")
    end <- start + attr(rr, "capture.length") - 1L
    substring(xx, start, end)
  }, xx = x, rr = rx, SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

.fastaComments2DataFrame <- function(x) {
  x <- .fastaCommentParser(x)
  DataFrame(DB = Rle(factor(x[, 1L])),
            AccessionNumber = x[, 2L],
            EntryName = x[, 3L],
            IsoformName = Rle(x[, 4L]),
            ProteinName = x[, 5L],
            OrganismName = Rle(factor(x[, 6L])),
            GeneName = x[, 7L],
            ## Levels of evidence:
            ## http://www.uniprot.org/manual/protein_existence
            ProteinExistence = Rle(factor(x[, 8L],
                                      labels = c("Evidence at protein level",
                                                 "Evidence at transcript level",
                                                 "Inferred from homology",
                                                 "Predicted",
                                                 "Uncertain"),
                                      levels = 1L:5L)),
            SequenceVersion = Rle(x[, 9L]))
}

