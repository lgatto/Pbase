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
\alias{length,Proteins-method}
\alias{mcols,Proteins-method}
\alias{metadata,Proteins-method}
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
    \item{\code{mcols}:}{\code{\linkS4class{DataFrame}}, containing metadata for
        each protein per row, e.g. \code{filenames}.}
    \item{\code{pfeatures}:}{\code{\linkS4class{IRanges}}.}
    \item{\code{seq}:}{\code{\linkS4class{AAStringSet}}, containing the sequence
        information of each protein.}
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

    \item{metadata}{\code{signature(x = "Proteins")}:
      Returns a \code{list} of global metadata.}

    \item{mcols}{\code{signature(x = "Proteins")}:
      Returns a \code{\linkS4class{DataFrame}} of metadata.}

    \item{length}{\code{signature(x = "Proteins")}:
      Returns the number of proteins.}

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
