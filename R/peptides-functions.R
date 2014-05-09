#' calculates the monoisotopic mass for peptide sequences
#' @param a character vector
#' @return a named doubled vector
#' @noRd
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

#' calculates the isoelectric point
#' @param a character vector
#' @return a named doubled vector
#' @references
#' Moore, Dexter S.
#' "Amino acid and peptide net charges: a simple calculational procedure."
#' Biochemical Education 13.1 (1985): 10-11.
#' http://dx.doi.org/10.1016/0307-4412(85)90114-1
#' @noRd
#.calculateIsoelectricPoint <- function(x) {
#}
#

#' test peptides for some properties
#' @param x character, AAString, AAStringSet: sequence
#' @param mass double, length == 2, mass range [Da]
#' @param length double, length == 2, length range
#' @return logical vector, TRUE if the peptides fulfills all criteria otherwise
#' FALSE
#' @noRd
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

    valid
}

