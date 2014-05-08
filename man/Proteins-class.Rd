\name{Proteins-class}
\Rdversion{1.1}
\docType{class}

\alias{Proteins-class}
\alias{class:Proteins}
\alias{Proteins}

\alias{[,Proteins,ANY,ANY-method}
\alias{[[,Proteins,ANY,ANY-method}
\alias{aa,Proteins-method}
\alias{aa}
\alias{ametadata,Proteins-method}
\alias{ametadata}
\alias{acols}
\alias{metadata,Proteins-method}
\alias{pfeatures,Proteins-method}
\alias{pfeatures}
\alias{plot,Proteins,missing-method}
\alias{pmetadata,Proteins-method}
\alias{pmetadata}
\alias{pcols}
\alias{pranges,Proteins-method}
\alias{pranges}

\alias{length,Proteins-method}
\alias{seqnames,Proteins-method}
\alias{show,Proteins-method}

\alias{proteinCoverage,Proteins,mzID-method}
\alias{proteinCoverage}
\alias{cleave,Proteins-method}

\title{The \code{Proteins} Class for Proteomics Data And Meta-Data}

\description{
  The \code{Proteins} class encapsulates data and meta-data for proteomics
  experiments, as described in the \code{slots} section.
}

\section{Objects from the Class}{

  Objects can be created using its constructor \code{Proteins}. The
  constructor either takes a \code{fasta} file name as first argument
  or, alternatively a named \code{uniprotIda} argument with valid
  UniProt accession numbers.
  
}

\section{Slots}{
  \describe{
    \item{\code{metadata}:}{Object of class \code{"list"} containing
      global metadata.
    }
    \item{\code{aa}:}{Object of class \code{"AAStringSet"} storing the
      protein sequences.
    }
    \item{\code{pranges}:}{Object of class
      \code{"CompressedIRangesList"} containg protein feature ranges
      such as theoretical (obtain by in silico cleavage) or observed
      peptides.
    }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      documenting the class verions. Intended for developer use and
      debugging.
    }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  \describe{


    \item{aa}{\code{signature(x = "Proteins")}: Returns an
      \code{\linkS4class{AAStringSet}} instance representing the
      sequences of the proteins.}
    
    \item{pfeatures}{\code{signature(x = "Proteins")}: ... }    
    \item{pranges}{\code{signature(x = "Proteins")}: ... }

    \item{metadata}{\code{signature(x = "Proteins")}: Returns a
      \code{list} of global metadata of the instance \code{x}, including
      data of instance creation or, if created from a set of UnitProt
      identifiers (see constructors above), the UniProt version and
      \code{UnitProt.WS} version number. }
    
    \item{ametadata}{\code{signature(x = "Proteins")}: Returns a
      \code{\linkS4class{DataFrame}} of protein metadata. }    
    \item{acols}{\code{signature(x = "Proteins")}: See
      \code{ametadata}. }    
    \item{pmetadata}{\code{signature(x = "Proteins")}: Returns a
      \code{list} of feature metadata. }    
    \item{pcols}{\code{signature(x = "Proteins")}: See
      \code{pmetadata}. } 
    
    \item{seqnames}{\code{signature(x = "Proteins")}: Returns the
      protein sequence names defined as UniProt accession numbers. }

    \item{length}{\code{signature(x = "Proteins")}: Returns the number
      of proteins. }
    
    \item{[}{\code{signature(x = "Proteins", i = "ANY", j = "missing")}:
      Creates a subset of the \code{Proteins} insance.}
      
    \item{[[}{\code{signature(x = "Proteins", i = "ANY", j = "missing")}:
      Returns an \code{\linkS4class{AAString}} instance representing the
      sequence of the selected protein.}

    \item{proteinCoverage}{\code{signature(x = "Proteins", y = "mzID")}:
      ... }    
    
    \item{cleave}{\code{signature(x = "Proteins")}: ... }
    
    \item{plot}{\code{signature(x = "Proteins", y = "missing")}: Plots
      all proteins and associated peptides using the
      \code{Gviz}/\code{Pviz} infrastructure. }

    \item{show}{\code{signature(object = "Proteins")}: Displays object
      summary as text. }
      
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

plot(p[1:5], from = 1, to = 30)
}


\keyword{classes}
