Pbase
=====

Manipulating and exploring protein and proteomics data.

Not sure if `Proteins` or `Pbase`; the latter fits well with `Pviz`.

## Ideas

### Protein coverage

`proteinCoverage` that return a summary (S4 class) of a protein
coverage given an experimental set of observed peptides (via an
`MSnExp`, `MSnSet` with an `fcol` or an `mzId` object) or a
theoretical set of peptides after in silico digestion (see also next
point). The return value would have its own `plot` method using
`Pviz`.

```r
setGeneric("proteinCoverage", function(x, y, ...), standardGeneric("proteinCoverage"))
setMethod("proteinCoverage", c("AAString", "MSnSet"), function(x, y, fcol = "pepseq") { ... } )
setMethod("proteinCoverage", c("AAString", "MSnExp"), function(x, y, fcol = "pepseq") { ... } )
setMethod("proteinCoverage", c("AAString", "mzID"), function(x, y) { ... } )
setMethod("proteinCoverage", c("AAString", "GRanges"), function(x, y) { ... } )
```
The above could also be vectorised to work with and `AAStringSet`.

### Assessing the redundancy of a protein fasta database

Given a protein fasta file, what is the maximal sensitivity that can
be expected from a mass spectrometry experiment with 0, 1,
... miscleavages. This should probably also include a filtering step
for peptide *flyability*.

### Selection of optimal heavy peptides for absolute quantitation

See Pavel's [idea](https://github.com/sgibb/cleaver/issues/5).

### Creating a link between arbitrary protein sequences and a genomic reference

If we want t o migrate towards a `GRanges` infrastructure, it is
important to be able to link proteins and full proteomics back to
genomic coordinates. I don't know if this mapping is provided for the
UniProt reference proteomes. 

See also the `sapFinder` Bioconductor packages.

### Protein domains

Maybe support for the annotation of detection of protein domains.

## Dependencies

The package should allow to easily interact with `AAString` and
`AAStringSet` instances, protein databases such as UniProt (and
possibly biomaRt) using protein identifiers, protein identification
results (`mzID` package and, later `mzR`) and possibly also `MSnExp`
and `MSnSet` instances. 

