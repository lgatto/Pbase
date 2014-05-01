Pbase
=====

Manipulating and exploring protein and proteomics data.

Uses `Pviz` for visualisation.

## Classes

```r
setClass("Protein",
         slots = list(
             fileName = "character", ## optional
             seq = "AAString",
             pfeatures = "IRanges",
             created = "character"
             UnitProtVersion = "character"),
         contains = "Versioned") ## also record Uniprot.WS and IRanges classes

setClass("Proteins",
         slots = list(
             fileName = "character", ## can be multiple files, optional
             seq = "AAStringSet",
             pfeatures = "IRanges",
             created = "character"
             UnitProtVersion = "character"),
         contains = "Versioned") ## also record Uniprot.WS and IRanges classes
```

## Constructor

`Protein("fastafile")` returns a `Protein` or `Proteins` instance, depending on the number of sequences.

`Protein("ids")` returns a `Protein` or `Proteins` instance, depending on the number of identifers.


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
## uses internal IRanges or, if absent, first cleaves the protein 
setMethod("proteinCoverage", c("Protein", "missing"), function(x, y) { ... } ) 
setMethod("proteinCoverage", c("Protein", "MSnSet"), function(x, y, fcol = "pepseq") { ... } )
setMethod("proteinCoverage", c("Protein", "MSnExp"), function(x, y, fcol = "pepseq") { ... } )
setMethod("proteinCoverage", c("Protein", "mzID"), function(x, y) { ... } )
```

And the same pf `Proteins`.

### Assessing the redundancy of a protein fasta database

Given a protein fasta file, what is the maximal sensitivity that can
be expected from a mass spectrometry experiment with 0, 1,
... miscleavages. This should probably also include a filtering step
for peptide *flyability*.

### Selection of obtimal heavy peptides for absolute quantitation

See Pavel's idea.

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

