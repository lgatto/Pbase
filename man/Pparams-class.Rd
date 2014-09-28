\name{Pparams-class}
\Rdversion{1.1}
\docType{class}
\alias{Pparams-class}
\alias{Pparams}
\alias{show,Pparams-method}



\title{Class \code{"Pparams"}}

\description{
  \code{Pbase} parametrisation infrastructure.
}

\section{Objects from the Class}{
  New \code{Pbase} parameters can be generated with the \code{Pparams()}
  constructor. \code{Pparams} instances control various aspects of
  \code{Pbase} functions, as described in the \emph{Slots} section
  below. If no parameters are passed to the respective functions,
  default values from \code{Pparams()} are used.
}


\section{Slots}{
  \describe{

    \item{\code{DbFormat}:}{The format of the protein sequence fasta
      database used to generate the \code{Proteins} object. Currently
      only \code{"UniProt"} is supported. \code{"RefSeq"} will be added
      as well as a mechanism to support arbitrary and custom fasta
      header.
    }

    \item{\code{IdFormat}:}{The format of the identification data files
      used to add \code{pfeatures} to \code{Protein}
      instances. Currently, \code{mzIdentML} is supported. }

    \item{\code{IdReader}:}{Package to be used to load the
      identification data. Currently one of \code{mzR} (via the
      \code{openIDfile} and \code{psms} functions) or \code{mzID} (via
      the \code{mzID} and \code{flatten} functions). Differences between
      these two architectures include the metadata available in the
      \code{Proteins}' \code{pfeatures}, speed and stability (\code{mzR}
      is much faster but less mature and currently susceptible to
      crashes). }

    \item{\code{verbose}:}{A \code{logical} defining if the various
      functions display messages (default) or remain silent. }
    
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "Pparams")}: ... }
	 }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
Pparams()
Pparams(IdReader = "mzID")

try(Pparams(IdReader = "mzid"))
}

\keyword{classes}
