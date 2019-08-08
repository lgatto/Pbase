\name{proteinCoding-methods}
\docType{methods}

\alias{proteinCoding}
\alias{proteinCoding-methods}
\alias{proteinCoding,GRanges-method}
\alias{proteinCoding,GRangesList-method}


\title{Only keep protein coding ranges}


\description{

  Removed all the ranges that are not protein coding. Typically used on
  the output of \code{\link{etrid2grl}}.

}


\section{Methods}{
  \describe{

    \item{\code{signature(object="GRanges", mcol="character",
	coding="character")}}{ Removes all the ranges that are not
	annotated as protein coding ranges, i.e. ranges whose
	\code{mcols()$mcol} is different from \code{coding}. The default
	values are \code{mcols()feature} and
	\code{"protein_coding"}. The method return the \code{GRanges}
	trimmed from all non-matching ranges.
    }

    \item{\code{signature(object="GRangesList", mcol="character",
	coding="character")}}{As above but for \code{GRanges} in a
      \code{GrangesList}.
    }
  }
}

\keyword{methods}
