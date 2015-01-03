#' Function to create heavy labeled peptides,
#' see https://github.com/sgibb/cleaver/issues/5 for details.
#' TODO: There should be a function to find the best labels for a given protein
#' automatically.
#' @param peptides a named character vector of peptides (names must correspond
#' to the protein accession numbers)
#' @param proteins a cleaved Proteins object
#' @param ... further arguments passed to .addOverhangs
#' @noRd
#'
.calculateHeavyLabels <- function(peptides, proteins, ...) {
    stopifnot(is(peptides, "character"))
    stopifnot(is(proteins, "Proteins"))

    if (is.null(names(peptides))) {
        stop("No names for ", sQuote("peptides"), " available!")
    }

    if (!isCleaved(proteins)) {
        stop("You have to ", sQuote("cleave"), " your proteins first!")
    }

    aa <- aa(proteins)
    pr <- pranges(proteins)
    pos <- .peptidePosition(peptides, aa)

    l <- vector(mode = "list", length = length(pos))

    ## TODO: Because .addOverhangs is only vectorized on peptide level (multiple
    ## peptides per protein) and not on protein level we need to loop over all
    ## proteins.
    for (i in seq_along(pos)) {
        pindex <- match(start(pos[[i]]), start(pr[[i]]))
        pindex <- pindex[!is.na(pindex)]

        if (length(pindex)) {
            l[[i]] <- cbind(Protein = names(aa)[i],
                            .addOverhangs(sequence = unlist(aa[[i]]),
                                          pindex = pindex,
                                          pranges = pr[[i]], ...),
                            stringsAsFactors = FALSE)
        }
    }

    d <- do.call(rbind, l)

    ## create labeled peptide string
    d$spikeTide <- apply(d, 1, function(x) {
        paste0(na.omit(c(x["N_overhang"], x["Peptide"], x["C_overhang"])),
               collapse=".")
    })

    return(d)
}

#' This functions creates peptide sequences for heavy labeled peptides. Sadly
#' this function is not completely vectorized yet.
#' TODO: implement method != "both"
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
.addOverhangs <- function(sequence, pindex, pranges, maxN = 20L,
                          nN = 4L, nC = 3L, endsWith = c("K", "R", "G"),
                          method = c("both", "N", "C")) {
    .isTRUE <- function(x) {
        sapply(x, isTRUE)
    }

    sequence <- as.character(sequence)
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

    iN <- which(w >= nN[1L])
    iC <- which(w >= nC[1L])

    ## handle too long sequences
    isTooLarge <- pw > (maxN - nN - nC)
    nN[isTooLarge] <- NA
    nC[isTooLarge] <- NA
    solution[isTooLarge] <- "no_overhangs"

    ## find first valid cleavage product
    ## shift iN left to avoid overflows
    before <- c(iN[1L], iN)[findInterval(pindex, iN, rightmost.closed = FALSE,
                                         all.inside = FALSE)]

    ## shift iC right to avoid overflows:
    after <- c(iC, tail(iC, 1L))[findInterval(pindex, iC,
                                              rightmost.closed = TRUE,
                                              all.inside = FALSE) + 1L]

    ## extend N- and C-terminus
    s <- pmax(end(pranges)[before] - (nN - 1L), 1L)
    e <- pmin(start(pranges)[after] + (nC - 1L), tail(end(pranges), 1L))

    ## test correct ending
    if (!is.null(endsWith)) {
        isCorrectEnding <- substring(sequence, e, e) %in% endsWith
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
    shortN <- .isTRUE(!isValidLength & cLength <= nC)
    shortC <- .isTRUE(!isValidLength & nLength <= nN)
    shortBoth <- .isTRUE(!isValidLength & !shortN & !shortC)

    s[shortN] <- (ps - (maxN - pw - nC))[shortN]
    e[shortC] <- (pe + (maxN - pw - nN))[shortC]
    s[shortBoth] <- (ps - nN)[shortBoth]
    e[shortBoth] <- (pe + nC)[shortBoth]
    solution[shortN] <- "N_overhang_shortened"
    solution[shortC] <- "C_overhang_shortened"
    solution[shortBoth] <- "both_overhangs_shortened"

    return(data.frame(Peptide = substring(sequence, ps, pe),
                      N_overhang = substring(sequence, s, ps - 1L),
                      C_overhang = substring(sequence, pe + 1L, e),
                      spikeTideResult = solution, stringsAsFactors = FALSE))
}

