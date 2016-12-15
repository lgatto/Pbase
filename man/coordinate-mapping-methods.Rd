\name{mapToGenome-methods}
\docType{methods}

\alias{mapToGenome}
\alias{mapToGenome-methods}
\alias{mapToGenome,ANY,ANY-method}
\alias{mapToGenome,Proteins,GRangesList-method}
\alias{mapToGenome,Proteins,EnsDb-method}

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
\S4method{mapToGenome}{Proteins,GRangesList}(x, genome, drop.empty.ranges = TRUE, ...)
\S4method{pmapToGenome}{Proteins,GRangesList}(x, genome,
drop.empty.ranges = TRUE, ...)
\S4method{mapToGenome}{Proteins,EnsDb}(x, genome, id = "name", idType =
"protein_id", drop.empty.ranges = TRUE, ...)
}


\arguments{

  \item{x}{\code{\linkS4class{Proteins}} object containing peptides
    \code{pranges} to be mapped. }

  \item{genome}{A \code{\linkS4class{GRangesList}} object used to map
    between \code{x} and the result. The ranges are typically created by
    the \code{\link{etrid2grl}} function. Alternatively, an
    \code{\link[ensembldb]{EnsDb}} object providing the required
    annotation for the mapping, i.e. the annotation of proteins to
    transcripts and the genomic coordinates of the transcripts' exons.}

  \item{drop.empty.ranges}{\code{TRUE} (default) or \code{FALSE}. Should
    empty ranges be dropped?}

  \item{id}{character(1) indicating which metadata columns in \code{x}
  provide the (protein) IDs for the mapping of proteins to
  transcripts. Can be the name of any columns in \code{acols(x)} or
  \code{"name"} in which \code{seqnames(x)} will be used.}

  \item{idType}{character(1) specifying the type of the IDs found in
  \code{id}. Supported are \code{"protein_id"} (the Ensembl protein ID),
  \code{"tx_id"} (Ensembl transcript ID) or \code{"uniprot_id"} (Uniprot
  ID).}

  \item{\ldots}{Additional parameters passed to inner
    functions. Currently ignored.  }
}

\details{

  \itemize{

    \item \code{mapToGenome} maps the \code{pranges(x)} to the ranges of
    \code{genome}. Unless \code{x} and \code{genome} are of length 1,
    both must be named and items of \code{x} are matched to items of
    \code{genome} using their respective names. Names that do not
    co-occur in \code{x} and \code{genome} are ignored. If we have

    \code{seqnames(x)}: \code{"A"}, \code{"B"} and \code{"C"}

    and

    \code{names(genome)}: \code{"C"}, \code{"A"}, \code{"a"},
    \code{"z"}, \code{"A"} and \code{"A"}.

    the names of the output will be

    \code{"A"}, \code{"A"}, \code{"A"} and \code{"C"}.

    The output is ordered by (1) \code{seqnames(x)} and (2) the order of
    the elements in \code{genome}.

    In case less than \code{length(x)} are mapped, as for \code{p["B"]}
    above, a message informs the user.

    \item \code{mapToGenome,Proteins,EnsDb} maps each of the
    \code{pranges(x)} ranges within the protein sequence to the
    corresponding genomic coordinates using annotations provided by the
    \code{\link[ensembldb]{EnsDb}} object. To enable the mapping the
    \code{Proteins} object has to provide IDs that can be used to
    identify the encoding transcript. Such IDs can be the Ensembl
    protein ID, the Uniprot ID or the Ensembl transcript ID. If a
    protein is annotated to multiple transcripts, the function selects
    the transcript which CDS length best matches the protein sequence
    length.

    The \code{mapToGenome,Proteins,EnsDb} method maps \code{pranges} of
    all proteins in the \code{Proteins} object to the genome. See
    examples below for more details.

    \item \code{pmapToGenome} is the element-wise (aka 'parallel')
    version of \code{mapToGenome}. The i-th \code{pranges(x)} is mapped
    to the i-th range in \code{genome}. \code{x} and \code{genome} must
    have the same length and do not need to be named (names are
    ignored).

  }

}

\seealso{
  \itemize{

    \item See \code{?\link[GenomicAlignments]{mapToAlignments}} in the
      \pkg{GenomicAlignments} package for mapping coordinates between
      reads (local) and genome (reference) space using a CIGAR
      alignment.

    \item See \code{?\link[GenomicFeatures]{mapToTranscripts}} in the
    \pkg{GenomicRanges} package for mapping coordinates between features
    in the transcriptome and genome space.

    \item The \code{\link{proteinCoding}} function to remove non-protein
    coding ranges before mapping peptides to their genomic coordinates.

    \item The \code{mapping} vignette for examples and visualisations.
  }

  See \code{\link{plotAsAnnotationTrack}} and
  \code{\link{plotAsAnnotationTrack}} for more details about the two
  plotting functions.

}



\value{

  A named \code{GRangesList} object, with names matching
  \code{names(genome)}. For \code{pmapToGenome}, the return value will
  have the same length as the inputs.

}

\author{Laurent Gatto, Johannes Rainer}

\examples{

data(p)

grl <- etrid2grl(acols(p)$ENST)
pcgrl <- proteinCoding(grl)

plotAsGeneRegionTrack(grl[[1]],
                      pcgrl[[1]])

mp <- mapToGenome(p[4], pcgrl[4])

plotAsAnnotationTrack(mp[[1]], pcgrl[[4]])

pmapToGenome(p, pcgrl)

#######
## mapToGenome,Proteins,EnsDb
## load an EnsDb object providing the required annotations
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## Map the pranges of all proteins in p to the genome providing the proteins'
## Uniprot IDs (being the 'names' of the Proteins object) for the mapping.
mp <- mapToGenome(p, edb, id = "name", idType = "uniprot_id")

}

\keyword{methods}
\keyword{utilities}
