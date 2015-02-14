\name{mapToGenome-methods}
\docType{methods}

\alias{mapToGenome}
\alias{mapToGenome-methods}
\alias{mapToGenome,ANY,ANY-method}
\alias{mapToGenome,Proteins,GenomicRanges-method}

\alias{pmapToGenome}
\alias{pmapToGenome-methods}
\alias{pmapToGenome,ANY,ANY-method}
\alias{pmapToGenome,Proteins,GRangesList-method}


\title{Map range coordinates between proteins and genome space}


\description{
  Map range coordinates between peptide features along proteins and
  genome (reference) space.
}


\usage{
\S4method{mapToGenome}{Proteins,GenomicRanges}(x, genome, ...)
\S4method{pmapToAlignments}{Proteins,GRangesList}(x, genome, ...)
}


\arguments{

  \item{x}{\code{\linkS4class{Proteins}} object contain peptides
    \code{pranges} to be mapped. Of length 1 of \code{genome} is a
    \code{\linkS4class{GenomicRanges}}. }

  \item{genome}{A named \code{\linkS4class{GenomicRanges}} or
    \code{\linkS4class{GRangesList}} object used to map between \code{x}
    and the result. The ranges are typically created by the
    \code{\link{etrid2grl}} function. }

  \item{\ldots}{Additional parameters passed to inner
    functions. Currently ignores.  }
}

\details{

  \code{mapToGenome} maps the \code{pranges(x)} to the ranges of
  \code{genome}. \code{x} is expected to be of length 1; any additional
  items will be ignored with a warning.

  \code{pmapToGenome} is the element-wise (aka 'parallel') version of
  \code{mapToGenome}. The i-th \code{pranges(x)} is mapped to the i-th
  range in \code{genome}. \code{x} and \code{genome} must have the same
  length. 

}

\seealso{
  \itemize{ 

    \item See \code{?\link[GenomicAlignments]{mapToAlignments}} in the
      \pkg{GenomicAlignments} package for mapping coordinates between
      reads (local) and genome (reference) space using a CIGAR
      alignment.

      \item See \code{?\link[GenomicFeatures]{mapToTranscripts}} in
      the \pkg{GenomicRanges} package for mapping coordinates between
      features in the transcriptome and genome space.

      \item The \code{\link{proteinCoding}} function to remove
      non-protein coding ranges before mapping peptides to their genomic
      coordinates.
    
  }
}


\value{

  An object of the same class as \code{genome}.

  When \code{genome} is a named \code{GrangesList}, then the retruned
  instance will be given the same names. 

}

\author{Laurent Gatto}

\examples{

data(p)

grl <- etrid2grl(acols(p)$ENST)
pcgrl <- proteinCoding(grl)
    
mapToGenome(p[4], pcgrl[[4]])

pmapToGenome(p, pcgrl)

}

\keyword{methods}
\keyword{utilities}
