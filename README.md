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

`Protein("ids")` returns a `Protein` or `Proteins` instance, depending on the number of identifiers.


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

### Mapping a Protein Sequence to a Genome Sequence

#### ACT3_DROME

Manual example for two UniProt - Identifiers:
1. [ACT3_DROME](http://www.uniprot.org/uniprot/P53501).
2. Follow the [UniParc link](http://www.uniprot.org/uniparc/UPI0000000EDE) in section *Sequences*.
3. Choose one EMBL entry, e.g. [ACV91653](http://www.ebi.ac.uk/ena/data/view/ACV91653).
4. Look for the Identifier in the *Sequence* section, e.g. *ACV91653.1* (fasta file comment, 3rd column).
5. Put these identifier into the search box of http://ensembl.org and choose the *Fruitfly* species.
6. Choose the first search result, e.g. [Act57B](http://www.ensembl.org/Drosophila_melanogaster/Gene/Summary?g=FBgn0000044&db=core).
7. Select the protein with the same length as the UniProt protein, e.g. [Act57B-RA/FBpp0071448](http://www.ensembl.org/Drosophila_melanogaster/Transcript/ProteinSummary?db=core;g=FBgn0000044;r=2R:16831533-16833945;t=FBtr0071519).
8. Choose the Sequence/Protein link in the left navigation bar to get the aminoacid sequence.
9. Compare both sequences.

Uniprot
```
>sp|P53501|ACT3_DROME Actin-57B OS=Drosophila melanogaster GN=Act57B PE=1 SV=1
MCDDEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
SKRGILTLKYPIEHGIITNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKM
TQIMFETFNSPAMYVAIQAVLSLYASGRTTGIVLDSGDGVSHTVPIYEGYALPHAILRLD
LAGRDLTDYLMKILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAAASTSLEKS
YELPDGQVITIGNERFRCPESLFQPSFLGMESCGIHETVYNSIMKCDVDIRKDLYANIVM
SGGTTMYPGIADRMQKEITSLAPSTIKIKIIAPPERKYSVWIGGSILASLSTFQQMWISK
EEYDESGPGIVHRKCF
```

Ensembl
```
MCDDEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
SKRGILTLKYPIEHGIITNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKM
TQIMFETFNSPAMYVAIQAVLSLYASGRTTGIVLDSGDGVSHTVPIYEGYALPHAILRLD
LAGRDLTDYLMKILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAAASTSLEKS
YELPDGQVITIGNERFRCPESLFQPSFLGMESCGIHETVYNSIMKCDVDIRKDLYANIVM
SGGTTMYPGIADRMQKEITSLAPSTIKIKIIAPPERKYSVWIGGSILASLSTFQQMWISK
EEYDESGPGIVHRKCF
```

#### ALBU_HUMAN

1. [ALBU_HUMAN](http://www.uniprot.org/uniprot/P02768).
2. Follow the [UniParc link](http://www.uniprot.org/uniparc/?query=uniprot:P02768&direct=yes) in section *Sequences*.
3. Choose one EMBL entry, e.g. CBL66721 is not working, choose the second one [ABJ16448](http://www.ebi.ac.uk/ena/data/view/ABJ16448).
4. Look for the Identifier in the *Sequence* section, e.g. *ABJ16448.1* (fasta file comment, 3rd column).
5. Put these identifier into the search box of http://ensembl.org and choose the *Human* species.
6. Choose the first search result, e.g. [ALB-001](http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000295897&db=core).
7. Select the protein with the same length as the UniProt protein, e.g. [ALB-001/ENST00000295897](http://www.ensembl.org/Homo_sapiens/Transcript/ProteinSummary?db=core;g=ENSG00000163631;r=4:74269956-74287129;t=ENST00000295897).
8. Choose the Sequence/Protein link in the left navigation bar to get the aminoacid sequence.
9. Compare both sequences.

#### Short cut

1. Visit the Uniprot page of the protein of interest.
2. Look for section *Sequence databases* and choose *EMBL* and an identifier.
3. Use this identifier to search on http://ensembl.org .

## Dependencies

The package should allow to easily interact with `AAString` and
`AAStringSet` instances, protein databases such as UniProt (and
possibly biomaRt) using protein identifiers, protein identification
results (`mzID` package and, later `mzR`) and possibly also `MSnExp`
and `MSnSet` instances.

