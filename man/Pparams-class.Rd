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
    \item{\code{DbFormat}:}{Object of class \code{"character"} ~~ }
    \item{\code{IdFormat}:}{Object of class \code{"character"} ~~ }
    \item{\code{IdReader}:}{Object of class \code{"character"} ~~ }
    \item{\code{verbose}:}{Object of class \code{"logical"} ~~ }
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
