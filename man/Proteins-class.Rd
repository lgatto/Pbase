\name{Proteins-class}
\Rdversion{1.1}
\docType{class}

\alias{Proteins-class}
\alias{class:Proteins}
\alias{Proteins}

\alias{Proteins,character,missing-method}
\alias{Proteins,missing,character-method}

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
\alias{pvarLabels,Proteins-method}
\alias{pvarLabels}
\alias{avarLabels,Proteins-method}
\alias{avarLabels}

\alias{length,Proteins-method}
\alias{seqnames,Proteins-method}
\alias{show,Proteins-method}


\alias{proteinCoverage,Proteins-method}
\alias{proteinCoverage}
\alias{proteotypic,Proteins-method}
\alias{proteotypic}
\alias{cleave,Proteins-method}
\alias{addIdentificationData,Proteins-method}

\title{The \code{Proteins} Class for Proteomics Data And Meta-Data}

\description{

  The \code{Proteins} class encapsulates data and meta-data for
  proteomics experiments. The class stores the protein sequences as well
  as specific subsets of interest, typically peptides, as ranges. The
  \code{Proteins} instances, the sequence and peptide slots are
  described by their respective metadata attributes.
  
}

\section{Objects from the Class}{

  Objects can be created using its constructor \code{Proteins}. The
  constructor either takes a \code{fasta} file name as first argument
  or, alternatively a named \code{uniprotIda} argument with valid
  UniProt accession numbers.
  
}

\section{Details}{
  An instance of class \code{Proteins} is characterised by one of
  multiple protein sequences that can be accessed as \code{AASrringSet}
  with the \code{aa} accessor. Sequence-specific annotation, such as
  accession numbers, protein and gene names, ... is available with the
  \code{ametadata} or \code{acols} methods. General metadata such as the
  data of creation of the instance are stored as a \code{list} returned
  by the \code{metadata} accessor.

  Each sequence of a \code{Proteins} instance can also be characterised
  by a set of specific ranges describing peptides of interest. These
  \emph{peptide features} can be extracted as an \code{AAStringSetList},
  where each sequence features 0 or more peptide features. These
  peptides features are encode as ranges along the original proteins
  sequences (a \code{list} of \code{IRanges}) that can be extracted with
  the \code{pranges} function. These peptide features have their own
  metadata describing for example peptide identification scores, number
  of missed cleavages, ... available with the \code{pmetadata} or
  \code{pcols} methods.
  


}


\section{Slots}{
  \describe{
    \item{\code{metadata}:}{ Object of class \code{"list"} containing
      global metadata, accessed with \code{metadata}.
    }
    \item{\code{aa}:}{ Object of class \code{"AAStringSet"} storing the
      protein sequences, accessed with \code{aa}.
    }
    \item{\code{pranges}:}{ Object of class
      \code{"CompressedIRangesList"} containg protein feature ranges
      such as theoretical (obtain by in silico cleavage) or observed
      peptides. Accessed as an \code{\linkS4class{IRangesList}} with
      \code{pranges} or and \code{\linkS4class{AAStringSetList}} with
      \code{pfeatres}. 
    }
    \item{\code{.__classVersion__}:}{ Object of class \code{"Versions"}
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

    \item{avarLabels}{\code{signature(x = "Proteins")}: Returns the
      names of the sequences metadata. }
    \item{pvarLabels}{\code{signature(x = "Proteins")}: Returns the
      names of the peptide feature metadata. }
    
    \item{seqnames}{\code{signature(x = "Proteins")}: Returns the
      protein sequence names defined as UniProt accession numbers. }

    \item{length}{\code{signature(x = "Proteins")}: Returns the number
      of proteins. }
    
    \item{[}{\code{signature(x = "Proteins", i = "ANY", j = "missing")}:
      Creates a subset of the \code{Proteins} insance.}
      
    \item{[[}{\code{signature(x = "Proteins", i = "ANY", j = "missing")}:
      Returns an \code{\linkS4class{AAString}} instance representing the
      sequence of the selected protein.}

    \item{proteinCoverage}{\code{signature(x = "Proteins")}: ... }

    \item{proteotypic}{\code{signature(x = "Proteins")}: ... }
    
    \item{cleave}{\code{signature(x = "Proteins")}: ... }

    \item{addIdentificationData}{\code{signature(object = "Proteins")}: ... }
    
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

pp <- cleave(p[1:100])
pp <- proteotypic(pp)
pp
pcols(pp[1:2])
plot(pp[1:2], from = 20, to = 30)

pp <- proteinCoverage(pp)
avarLabels(pp)
pp
}


\keyword{classes}
