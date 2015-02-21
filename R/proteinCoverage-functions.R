proteinCoverage <- function(x) {
    stopifnot(is(x, "Proteins"))   
    pattern <- pranges(x)
    subject <- aaranges(x, unshift = TRUE)
    cvg <- .proteinCoverage(pattern, subject)    
    addacol(x, "Coverage", cvg)
}


### Calculates the coverage of patterns in subjects based on the width of the
### patterns. This method assumes that all patterns lie within a subject and the
### the names of patterns and subjects containing the correct accession.
### Caution: This is based purley on IRanges. No sequence based checks are
### involved! You have to make sure that you compare compareable sequences.
### @param pattern IRangesList of interest (from mzID, ...)
### @param subject IRangesList to compare against (mostly from Proteins)
### @return double, NA for patterns that are outside of subjects.
.proteinCoverage <- function(pattern, subject) {
    stopifnot(is(pattern, "IRangesList"))
    stopifnot(is(subject, "IRangesList"))

    if (is.null(names(pattern))) {
        stop("No names for ", sQuote("pattern"), " available!")
    } else if (anyDuplicated(names(pattern))) {
        stop("No duplicated names for ", sQuote("pattern"), " allowed!")
    }
    if (is.null(names(subject))) {
        stop("No names for ", sQuote("subject"), " available!")
    } else if (anyDuplicated(names(subject))) {
        stop("No duplicated names for ", sQuote("subject"), " allowed!")
    }

    ## combine overlapping ranges
    pattern <- reduce(pattern)

    ## reorder patterns by subject names (introduces NULL for non matches)
    orderedPattern <- as.list(pattern)[names(subject)]

    ## replace NULL by empty IRanges()
    orderedPattern[unlist(lapply(orderedPattern, is.null))] <- IRanges()
    orderedPattern <- IRangesList(orderedPattern)

    ## test that all patterns lie within the subject
    patternEnd <- as.integer(lapply(end(orderedPattern), tail, n = 1L))
    subjectEnd <- as.integer(lapply(end(subject), tail, n = 1L))
    isOutside <- which(patternEnd > subjectEnd)

    ## calculate coverage
    cvg <- as.double(sum(width(orderedPattern))/width(subject))
    cvg[isOutside] <- NA
    names(cvg) <- names(subject)

    return(cvg)
}

