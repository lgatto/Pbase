##' calculates the monoisotopic mass for peptide sequences
##' @param a character vector
##' @return a named doubled vector
##' @noRd
.calculateMolecularWeight <- function(x) {
    water <- sum(.get.atomic.mass()[c("H", "H", "O")])

    aamass <- setNames(.get.amino.acids()$ResidueMass,
                       .get.amino.acids()$AA)

    ppmass <- unlist(lapply(.singleAA(x), function(p) {
        sum(aamass[p])
    }))
    ppmass <- ppmass + water

    .setNames2(ppmass, x)
}

##' calculates the isoelectric point
##' @param a character vector
##' @return a named doubled vector
##' @references
##' Moore, Dexter S.
##' "Amino acid and peptide net charges: a simple calculational procedure."
##' Biochemical Education 13.1 (1985): 10-11.
##' http://dx.doi.org/10.1016/0307-4412(85)90114-1
##' @noRd
## .calculateIsoelectricPoint <- function(x) {
## }

##' test peptides for some properties
##' @param x character, AAString, AAStringSet: sequence
##' @param mass double, length == 2, mass range [Da]
##' @param len double, length == 2, length range
##' @return logical vector, TRUE if the peptides fulfills all criteria otherwise
##' FALSE
##' @noRd
.isValidPeptide <- function(x, mass = NULL, len = NULL
                            ## to be extended (isoelectric point, ...)
                            ) {
    ## initialize with TRUE
    valid <- !logical(length(x))

    if (!is.null(mass)) {
        valid <- valid & .isInRange(.calculateMolecularWeight(x), range = mass)
    }
    if (!is.null(len)) {
        valid <- valid & .isInRange(nchar(x), range = len)
    }

    .setNames2(valid, x)
}

##' calculates IRanges for peptides "pattern" in a protein "subject"
##' @param pattern named character, AAString, AAStringSet, AAStringSetList
##' @param subject named character, AAString, AAStringSet
##' @return a named IRangesList The index of the peptides
##' (PeptideIndex) and the proteins (ProteinIndex) is stored in
##' \code{IRangesList@@elementMetaData} as columns of a DataFrame.
##' @noRd
.peptidePosition <- function(pattern, subject) {
    if (is.null(names(pattern))) {
        stop("No names for ", sQuote("pattern"), " available!")
    }
    if (is.null(names(subject))) {
        stop("No names for ", sQuote("subject"), " available!")
    } else if (anyDuplicated(names(subject))) {
        stop("No duplicated names for ", sQuote("subject"), " allowed!")
    }

    if (is(pattern, "AAStringSetList")) {
        pattern <- unlist(pattern)
    }

    proteinIndex <- match(names(pattern), names(subject))

    l <- vector(mode = "list", length = length(pattern))

    ## convert subject from AAStringSet to character;
    ## This decreases the running time of the for loop dramatically (because it
    ## avoids a lot of standardGeneric dispatching calls and conversions to
    ## character).
    subject <- setNames(as.character(subject), names(subject))
    pattern <- as.character(pattern)

    for (i in seq(along = l)) {
        l[[i]] <- as.vector(gregexpr(pattern = pattern[i],
                                     text = subject[proteinIndex[i]],
                                     fixed = TRUE)[[1L]])
    }

    matches <- unlist(l)
    nmatches <- elementNROWS(l)
    nchars <- rep.int(nchar(pattern), nmatches)

    isMatch <- matches != -1

    matches <- matches[isMatch]
    nmatches <- nmatches[isMatch]
    nchars <- nchars[isMatch]

    ir <- IRanges(start = matches, width = nchars)

    mcols(ir) <-
        DataFrame(PeptideIndex = Rle(rep.int(seq_along(pattern)[isMatch],
                                             nmatches)),
                  ProteinIndex = Rle(rep.int(proteinIndex[isMatch], nmatches)))
    split(ir, f = names(subject)[as.integer(mcols(ir)$ProteinIndex)])
}

