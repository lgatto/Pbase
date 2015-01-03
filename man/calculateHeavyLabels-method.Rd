\name{calculateHeavyLabels-methods}
\docType{methods}
\alias{calculateHeavyLabels}
\alias{calculateHeavyLabels,Proteins-method}
\title{Calculate heavy labeled peptides.}
\description{
This method calculate heavy labeled peptides for proteins stored in a
\code{\linkS4class{Proteins}} object.
}

\usage{
\S4method{calculateHeavyLabels}{Proteins}(object, peptides,
                                          maxN = 20L, nN = 4L, nC = 3L,
                                          endsWith = c("K", "R", "G"))
}
\arguments{
  \item{object}{A \code{\linkS4class{Proteins}} object.}
  \item{peptides}{A named \code{character} vector containing the peptides of
      interest. The names must match the UniProt accession numbers of the
      proteins in \code{object}.}
  \item{maxN}{\code{interger}, maximal length of the heavy labeled peptide.}
  \item{nN}{\code{interger}, minimal number of amino acids at the N terminus.}
  \item{nC}{\code{interger}, minimal number of amino acids at the C terminus.}
  \item{endsWith}{A \code{character} vector containing the allowed amino acids
      at the end of the resulting sequence (every peptide that doesn't end with
      one of these amino acids has to be one amino acid shorter as
      \code{maxN}).}
}
\details{
The digestion efficiency with enzymes like trypsin is below 100\%. That's why
spiked-in peptides for labeled quantitation have to follow the same digestion
rules as the peptides of interest. Therefore it is necessary to extend the
peptides of interest by a few amino acids on the N- and C-terminus. These
extensions should not be a cleavage point of the used enzym. This methods
provides an easy interface to find the sequences for heavy labeled peptides that
could be used as spike-ins for the peptides of interest. Please see the
references for a more detailed discussion.
}
\value{
Returns a \code{data.frame} with 6 columns:
    \itemize{
        \item{Protein}{The Protein accession number.}
        \item{Peptide}{The peptide of interest.}
        \item{N_overhang}{The added sequence of the N-terminus.}
        \item{C_overhang}{The added sequence of the C-terminus.}
        \item{spikeTideResult}{A short description of the used creation rule.}
        \item{spikeTide}{The heavy labeled peptide that represents the peptide
            of interest best.}
    }
}
\references{
The complete description of the issue:
\url{https://github.com/sgibb/cleaver/issues/5}

Kito, Keiji, et al.
A synthetic protein approach toward accurate mass spectrometric quantification
of component stoichiometry of multiprotein complexes.
Journal of proteome research 6.2 (2007): 792-800.
\url{http://dx.doi.org/10.1021/pr060447s}
}
\author{
Sebastian Gibb \email{mail@sebastiangibb.de} and
Pavel Shliaha \email{pvs22@cam.ac.uk}
}
\seealso{
\code{\linkS4class{Proteins}}
}
\examples{
## example protein database
data(p, package = "Pbase")

## digest proteins into peptides
cleavedProteins <- cleave(p)

## find spike-ins for the peptides of interest
calculateHeavyLabels(cleavedProteins,
                     peptides = c(A4UGR9 = "MEGFHIK",
                                  A4UGR9 = "QGNMYTLSK",
                                  A6H8Y1 = "GSTASNPQR"))
}
\keyword{methods}
