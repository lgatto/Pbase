## GenBank GI number changes each time the sequence changes; the
## accession number remains the same.
## http://www.ncbi.nlm.nih.gov/Sitemap/sequenceIDs.html
.fastaCommentParser_refseq <- function(x) {
    ## >gi|262118207|ref|NM_000202.5| Homo sapiens iduronate 2-sulfatase (IDS), transcript variant 1, mRNA
    ## >gi|530357116|ref|XP_005275933.1| PREDICTED: peptidyl-prolyl cis-trans isomerase A-like [Homo sapiens]
    ## short perl explanation:
    ## (?<NAME>...) creates a named group
    ## (?:...) creates a group but doesn't list it in the results
    ## (...)? optional group
    ## \\s matches spaces

    rx <- gregexpr(pattern = paste0(
                       ## "gi"
                       "^>?(?:.+)\\|",
                       ## gi number, nums only
                       "(?<GI>[0-9]+)\\|",
                       ## ref (for RefSeq)
                       "(?:.+\\|)",
                       ## Accession number
                       "(?<AN>[A-Z_0-9-\\.]+)\\|\\s",
                       ## ProteinName, e.g. Lipoprotein
                       "(?<PN>.+)$"), text = x, perl = TRUE)
    x <- do.call(rbind, .perlRxSubstring(x, rx))
    x[!nchar(x)] <- NA_character_
    x
}

.isRefSeqAccessionNumber <- function(x) {
    grepl(paste0("^[NYX]P_[0-9]+[\\.]?[0-9]+?"), x)
}



