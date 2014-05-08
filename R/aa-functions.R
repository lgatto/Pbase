#' calculates the monoisotopic mass for peptide sequences
#' @param a character vector
#' @return a named doubled vector
.calculateMolecularWeight <- function(x) {
    water <- sum(get.atomic.mass()[c("H", "H", "O")])

    aamass <- setNames(get.amino.acids()$ResidueMass,
                       get.amino.acids()$AA)

    ppmass <- unlist(lapply(strsplit(as.character(x), ""), function(p) {
        sum(aamass[p])
    }))
    ppmass <- ppmass + water

    if (is.null(names(x))) {
        names(ppmass) <- as.character(x)
    } else {
        names(ppmass) <- names(x)
    }

    ppmass
}

