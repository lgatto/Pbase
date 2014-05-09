#' calculates the monoisotopic mass for peptide sequences
#' @param a character vector
#' @return a named doubled vector
#' @noRd
.calculateMolecularWeight <- function(x) {
    water <- sum(.get.atomic.mass()[c("H", "H", "O")])

    aamass <- setNames(.get.amino.acids()$ResidueMass,
                       .get.amino.acids()$AA)

    ppmass <- unlist(lapply(strsplit(as.character(x), ""), function(p) {
        sum(aamass[p])
    }))
    ppmass <- ppmass + water

    .setNames2(ppmass, x)
}

