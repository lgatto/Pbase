##' A function to calculate heavy labeled peptides for proteins stored
##' in a \code{\linkS4class{Proteins}} object.
##'
##' The digestion efficiency with enzymes like trypsin is below
##' 100\%. That's why spiked-in peptides for labeled quantitation have
##' to follow the same digestion rules as the peptides of
##' interest. Therefore it is necessary to extend the peptides of
##' interest by a few amino acids on the N- and C-terminus. These
##' extensions should not be a cleavage point of the used enzym. This
##' methods provides an easy interface to find the sequences for heavy
##' labeled peptides that could be used as spike-ins for the peptides
##' of interest. Please see the references for a more detailed
##' discussion.
##'
##' TODO: There should be a function to find the best labels for a
##' given protein automatically.
##' 
##' @title Calculate heavy labeled peptides
##' @param proteins A \code{\linkS4class{Proteins}} object.
##' @param peptides A named \code{character} vector containing the
##' peptides of interest. The names must match the UniProt accession
##' numbers of the proteins in \code{object}.
##' @param maxN An \code{integer}, maximal length of the heavy labeled
##' peptide.
##' @param nN An \code{integer}, minimal number of amino acids at the
##' N terminus.
##' @param nC An \code{integer}, minimal number of amino acids at the
##' C terminus.
##' @param endsWith A \code{character} vector containing the allowed
##' amino acids at the end of the resulting sequence (every peptide
##' that doesn't end with one of these amino acids has to be one amino
##' acid shorter as \code{maxN}).
##' @param ... Additional parameters passed to \code{.addOverhangs}.
##' @return A \code{data.frame} with 6 columns:
##' \itemize{
##'   \item{Protein}{The Protein accession number.}
##'   \item{Peptide}{The peptide of interest.}
##'   \item{N_overhang}{The added sequence of the N-terminus.}
##'   \item{C_overhang}{The added sequence of the C-terminus.}
##'   \item{spikeTideResult}{A short description of the used creation rule.}
##'   \item{spikeTide}{The heavy labeled peptide that represents the
##'         peptide of interest best.}
##' }
##' @author Sebastian Gibb <mail@@sebastiangibb.de> and Pavel Shliaha
##' @references
##' The complete description of the issue:
##' \url{https://github.com/sgibb/cleaver/issues/5}
##'
##' Kito, Keiji, et al. A synthetic protein approach toward accurate
##' mass spectrometric quantification of component stoichiometry of
##' multiprotein complexes. Journal of proteome research 6.2 (2007):
##' 792-800. \url{http://dx.doi.org/10.1021/pr060447s}
##' @examples
##' ## example protein database
##' data(p, package = "Pbase")
##'
##' ## digest proteins into peptides
##' cleavedProteins <- cleave(p)
##'
##' ## find spike-ins for the peptides of interest
##' calculateHeavyLabels(cleavedProteins,
##'                       peptides = c(A4UGR9 = "MEGFHIK",
##'                                    A4UGR9 = "QGNMYTLSK",
##'                                    A6H8Y1 = "GSTASNPQR"))
calculateHeavyLabels <-
              function(proteins, peptides,
                       maxN = 20L, nN = 4L, nC = 3L,
                       endsWith = c("K", "R", "G"),
                       ...) {
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

##' This functions creates peptide sequences for heavy labeled peptides. Sadly
##' this function is not completely vectorized yet.
##' TODO: implement method != "both"
##' See https://github.com/sgibb/cleaver/issues/5 for details.
##' @param sequence character/AAString/AAStringSet, protein sequence(s)
##' @param pindex index (position) of the peptide(s) in the cleaved protein(s)
##' @param pranges IRangesList, cleaved protein(s)
##' @param maxN integer, maximal length of the heavy labeled peptide
##' @param nN integer, minimal number of AA at the N terminus
##' @param nC integer, minimal number of AA at the C terminus
##' @param endsWith character, accepted ending AA (every peptide that doesn't end
##' with these AA has to be one AA shorter).
##' @param method character, should the N/C terminus prefered or should they
##' treated equaly; not implemented yet
##' @noRd
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
    if (length(before)) {
        s <- pmax(end(pranges)[before] - (nN - 1L), 1L)
    } else {
        s <- 1L
    }
    if (length(after)) {
        e <- pmin(start(pranges)[after] + (nC - 1L), tail(end(pranges), 1L))
    } else {
        e <- tail(end(pranges), 1L)
    }

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

