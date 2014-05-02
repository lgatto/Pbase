\name{Proteins-class}
\Rdversion{1.1}
\docType{class}
\concept{Proteins}

\alias{Proteins-class}
\alias{class:Proteins}
\alias{Proteins}

\alias{Proteins,character,missing-method}
\alias{Proteins,missing,character-method}
\alias{[,Proteins-method}
\alias{[[,Proteins-method}
\alias{aa}
\alias{aa,Proteins-method}
\alias{acols}
\alias{acols,Proteins-method}
\alias{ametadata}
\alias{ametadata,Proteins-method}
\alias{length,Proteins-method}
\alias{metadata,Proteins-method}
\alias{pcols}
\alias{pcols,Proteins-method}
\alias{pmetadata}
\alias{pmetadata,Proteins-method}
\alias{show,Proteins-method}

\title{The 'Proteins' Class for Proteomics Data And Meta-Data}

\description{
  The \code{Proteins} class encapsulates data and meta-data for proteomics
  experiments, as described in the \code{slots} section.
}

\section{Objects from the Class}{
  Objects can be created using its constructor \code{Proteins}.
}

\section{Slots}{
  \describe{
    \item{\code{metadata}:}{\code{list}, containing global metadata.}

    \item{\code{aa}:}{\code{\linkS4class{AAStringSet}}, containing the sequence
        information of each protein.}
    \item{\code{pfeatures}:}{\code{\linkS4class{IRangesList}}.}
  }
}

\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  \describe{
    \item{Proteins}{\code{signature(file = "character",
        uniprotIds = "missing")}:
      Creates an \code{\linkS4class{Proteins}} instance from one or multiple
      fasta files.}

    \item{Proteins}{\code{signature(file = "missing",
        uniprotIds = "character")}:
      \bold{NOT IMPLEMENTED YET}
    }

    \item{[}{\code{signature(x = "Proteins", ...)}:
      Creates a subset of \code{\linkS4class{Proteins}}.}

    \item{[[}{\code{signature(x = "Proteins", ...)}:
      Returns an \code{\linkS4class{AAString}} instance representing the
      sequence of the selected protein.}

    \item{aa}{\code{signature(x = "Proteins")}:
      Returns an \code{\linkS4class{AAStringSet}} instance representing the
      sequences of the proteins.}

    \item{acols}{\code{signature(x = "Proteins")}:
      See \code{\link{ametadata}}.}

    \item{ametadata}{\code{signature(x = "Proteins")}:
      Returns a \code{\linkS4class{DataFrame}} of protein metadata.}

    \item{metadata}{\code{signature(x = "Proteins")}:
      Returns a \code{list} of global metadata.}

    \item{length}{\code{signature(x = "Proteins")}:
      Returns the number of proteins.}

    \item{pcols}{\code{signature(x = "Proteins")}:
      See \code{\link{ametadata}}.}

    \item{pmetadata}{\code{signature(x = "Proteins")}:
      Returns a \code{\linkS4class{DataFrame}} of peptide metadata.}

    \item{show}{\code{signature(object = "Proteins")}:
      Displays object content as text. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk> and Sebastian Gibb <mail@sebastiangibb.de>
}

\examples{
fastaFiles <- list.files(system.file("extdata", package = "Pbase"),
                         pattern = "fasta", full.names = TRUE)
p <- Proteins(fastaFiles)
p
}

\keyword{classes}
