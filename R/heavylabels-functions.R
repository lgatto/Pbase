#' Function to create heavy labeled peptides,
#' see https://github.com/sgibb/cleaver/issues/5 for details.
#' TODO: There should be a function to find the best labels for a given protein
#' automatically.
#' @noRd
#'
dummyAddOverhang <- function(peptides, proteins, ...) {
    stopifnot(is(peptides, "character"))
    stopifnot(is(proteins, "Proteins"))

    if (is.null(names(peptides))) {
        stop("No names for ", sQuote("peptides"), " available!")
    }

    proteins <- proteins[unique(names(peptides))]

    pos <- .peptidePosition(peptides, aa(proteins))
    pr <- pranges(proteins)
    aa <- aa(proteins)

    shiftBy <- c(0L, cumsum(head(nchar(aa), -1L)))
    pos <- .flatIRangesList(pos, shift = TRUE, shiftBy = shiftBy)
    pr <- .flatIRangesList(pr, shift = TRUE, shiftBy = shiftBy)

    pindex <- match(start(pr), start(pos))
    pindex <- pindex[!is.na(pindex)]

    print(pindex)

    .addOverhang(sequence = unlist(aa), pindex = pindex, pranges = pr, ...)
}

#' This functions creates peptide sequences for heavy labeled peptides.
#' See https://github.com/sgibb/cleaver/issues/5 for details.
#' @param sequence character/AAString/AAStringSet, protein sequence(s)
#' @param pindex index (position) of the peptide(s) in the cleaved protein(s)
#' @param pranges IRangesList, cleaved protein(s)
#' @param maxN integer, maximal length of the heavy labeled peptide
#' @param nN integer, minimal number of AA at the N terminus
#' @param nC integer, minimal number of AA at the C terminus
#' @param endsWith character, accepted ending AA (every peptide that doesn't end
#' with these AA has to be one AA shorter).
#' @param method character, should the N/C terminus prefered or should they
#' treated equaly; not implemented yet
#' @noRd
.addOverhang <- function(sequence, pindex, pranges, maxN = 20L,
                         nN = 4L, nC = 3L, endsWith = c("K", "R", "G"),
                         method = c("both", "N", "C")) {
    stopifnot(is(pranges, "IRanges"))

    method <- match.arg(method)

    ## TODO:
    if (method != "both") {
        stop("method != \"both\" isn't supported yet")
    }

    maxN <- rep_len(maxN, length(pindex))
    nN <- rep_len(nN, length(pindex))
    nC <- rep_len(nC, length(pindex))
    solution <- character(length(pindex))

    w <- width(pranges)
    ps <- start(pranges)[pindex]
    pe <- end(pranges)[pindex]
    pw <- w[pindex]

    iN <- which(w >= nN)
    iC <- which(w >= nC)

    ## handle too long sequences
    isTooLarge <- pw > (maxN - nN - nC)
    nN[isTooLarge] <- NA
    nC[isTooLarge] <- NA
    solution[isTooLarge] <- "no_overhangs"

    ## find first valid cleavage product
    before <- iN[findInterval(pindex, iN, rightmost.closed = FALSE,
                              all.inside = FALSE) - 1L]
    after <- iC[findInterval(pindex, iC, rightmost.closed = TRUE,
                             all.inside = FALSE) + 1L]

    ## extend N- and C-terminus
    s <- pmax(end(pranges)[before] - (nN - 1L), 1L)
    e <- pmin(start(pranges)[after] + (nC - 1L), tail(end(pranges), 1L))

    ## test correct ending
    if (!is.null(endsWith)) {
        isCorrectEnding <- substr(sequence, e, e) %in% endsWith
        maxN[!isCorrectEnding] <- maxN[!isCorrectEnding] - 1L
        e[!isCorrectEnding] <- e[!isCorrectEnding] - 1L
    }

    ## 1. test length
    newLength <- (e + 1L) - s
    isValidLength <- newLength <= maxN
    solution[isValidLength] <- "fully_representative"

    ## 2. test N-/C-terminus length
    cLength <- e - pe
    nLength <- ps - s

    ## 3. shorten N/C or both
    shortN <- !isValidLength & cLength <= nC
    shortC <- !isValidLength & nLength <= nN
    shortBoth <- !isValidLength & !shortN & !shortC

    s[shortN] <- ps - (maxN - pw - nC)
    e[shortC] <- pe + (maxN - pw - nN)
    s[shortBoth] <- ps - nN
    e[shortBoth] <- pe + nC
    solution[shortN] <- "N_overhang_shortened"
    solution[shortC] <- "C_overhang_shortened"
    solution[shortBoth] <- "both_overhangs_shortened"

    ## TODO: return something more useful
    #return(list(pranges=IRanges(s, e), solution=solution))
    return(data.frame(Peptide = substr(sequence, ps, pe),
                      N_overhang = substr(sequence, s, ps - 1L),
                      C_overhang = substr(sequence, pe + 1L, e),
                      spikeTideResult = solution, stringsAsFactors = FALSE))
}

