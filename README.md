Pbase
=====

Manipulating and exploring protein and proteomics data.

Uses `Pviz` for visualisation.

## Classes

See `AllClasses.R`

Current decision is to avoid `pranges` from multiple origins. It is
probably easier to manage this situation using different `Proteins`
instances. We can then think about *comparing* `Proteins` instances
(that have `identical(p1@aa, p2@aa)` and easily return *common*
sequences.

## Accessors

- `pfeatures` returns an `AAStringSetList` of peptides (using `extractAt`)
- `pranges` returns the `pranges` slot (`CompressedIRangesList`) - was `pfeatures`.
- `aa` returns the protein sequences as an `AAStringSet`.

## metadata

- global: slot `metadata`, accessor `metadata`
  - `created` character
  - `UniProtVersion` character
  - `Uniprot.WS`, optional, only if constructed via ids.


- `aa` metadata:
  - `mcols(.@seq)` with accessor `acols` and `ametadata`
  - `metadata(.@seq)` - ignore

- `pfeatures` metadata:
  - `lapply(.@pfeatures, mcols)` with accessor `pcols` `pmetadata`
  - `metadata(.@pfeatures)` - ignore

## Constructor

`Proteins("fastafile")` returns a `Proteins` instance.

- Definition of the UniProt fasta comment format:
    http://www.uniprot.org/help/fasta-headers

`Proteins("ids")` returns a `Proteins` instance, depending on the
number of identifiers.

## Ideas

### Protein coverage

`proteinCoverage` that return a summary (S4 class) of a protein
coverage given an experimental set of observed peptides (via an
`MSnExp`, `MSnSet` with an `fcol` or an `mzId` object) or a
theoretical set of peptides after in silico digestion (see also next
point). The return value would have its own `plot` method using
`Pviz`.

```s
setGeneric("proteinCoverage", function(x, y, ...), standardGeneric("proteinCoverage"))
## uses internal IRanges or, if absent, first cleaves the protein
setMethod("proteinCoverage", c("Proteins", "missing"), function(x, y) { ... } )
## later
## setMethod("proteinCoverage", c("Proteins", "MSnSet"), function(x, y, fcol = "pepseq") { ... } )
## setMethod("proteinCoverage", c("Proteins", "MSnExp"), function(x, y, fcol = "pepseq") { ... } )
##setMethod("proteinCoverage", c("Proteins", "mzID"), function(x, y) { ... } )
```

### Assessing the redundancy of a protein fasta database

Given a protein fasta file, what is the maximal sensitivity that can
be expected from a mass spectrometry experiment with 0, 1,
... miscleavages. This should probably also include a filtering step
for peptide *flyability*.

#### Flyability/Detectability

Some literature about estimating detectability:

- LogR: [Liu, Hui, et al. "The Prediction of Peptide Detectability in MS Data Analysis Using Logistic Regression." Bioinformatics and Biomedical Engineering,(iCBBE) 2011 5th International Conference on. IEEE, 2011.](http://dx.doi.org/10.1109/icbbe.2011.5780167)
- SVM: [Webb-Robertson, Bobbie-Jo M., et al. "A support vector machine model for the prediction of proteotypic peptides for accurate mass and time proteomics." Bioinformatics 24.13 (2008): 1503-1509.](http://dx.doi.org/10.1093/bioinformatics/btn218)
- NN: [Sanders, William S., et al. "Prediction of peptides observable by mass spectrometry applied at the experimental set level." BMC bioinformatics 8.Suppl 7 (2007): S23.](http://dx.doi.org/10.1186/1471-2105-8-S7-S23)
- Gausian Mixed Discrimination [Mallick, Parag, et al. "Computational prediction of proteotypic peptides for quantitative proteomics." Nature biotechnology 25.1 (2007): 125-131.](http://dx.doi.org/10.1038/nbt1275)

##### Liu et al. 2011:

Requirements for in-silico created peptides: `missedCleavages = 0:2`, `length(peptides) >= 6`, `mass(peptides) < 6000` (Da)

Logistic Regression based on Hydrophobicity, Isoelectric point, length,
molecular weight, average hydrophobicity, average isoelectric point

##### Webb-Robertson et al. 2007:

Requirements for in-silico created peptides: `missedCleavages = 0:2`, `length(peptides) >= 6`, `mass(peptides) < 6000` (Da)

35 features: length, weidght, # of (non-)polar, # of (un)charged, # of pos./neg. charged residues, hydrophobicity (different models), polarity (different models), bulkiness, AA singlet counts

##### Sanders et al. 2007

Requirements for in-silico created peptides: `length(peptides) >= 6`

Features: Length, Charge, Isoelectric Point, Molecular Weight, Hydropathicity, Counts of each AA (20 Features), Percent composition of each AA (20 Features), Percent of polar, psoitive, negative, hydrophobic AA

take-home-message: a model of one species/dataset could not be transfered to another dataset (without dramatically decreasing the performance)

##### Mallick et al. 2007

~1000 Features.

Some of the most discriminating properties:
Total/Average net/positive charge, hydrophobic moment, isoelectric point, Histidine composition

take-home-message: The model of one species is comparable to another if the evolutionary
distance is small (e.g. yeast and human) but you can't compare different devices/datasets (e.g. MALDI vs ESI)

###### Simple Rules

Mass: `500:4500`

http://www.nature.com/nbt/journal/v25/n1/extref/nbt1275-S5.pdf
http://ieeexplore.ieee.org/ielx5/5779756/5779971/5780167/html/img/5780167-fig-1-large.gif

Length: `5:40`

http://www.nature.com/nbt/journal/v25/n1/extref/nbt1275-S6.pdf
http://ieeexplore.ieee.org/ielx5/5779756/5779971/5780167/html/img/5780167-fig-1-large.gif

95% of all peptides are of length `5:30`:
http://www.nature.com/nbt/journal/v25/n1/extref/nbt1275-S24.pdf

Average Isoelectric point: `seq(0, 1.4)`
http://ieeexplore.ieee.org/ielx5/5779756/5779971/5780167/html/img/5780167-fig-1-large.gif

### Hydropathy/Hydrophobicity

http://web.expasy.org/tools/protparam/protparam-doc.html
http://web.expasy.org/compute_pi/pi_tool-doc.html
[Kyte, Jack, and Russell F. Doolittle. "A simple method for displaying the hydropathic character of a protein." Journal of molecular biology 157.1 (1982): 105-132.](http://dx.doi.org/10.1016/0022-2836(82)90515-0)

### Selection of optimal heavy peptides for absolute quantitation

See Pavel's [idea](https://github.com/sgibb/cleaver/issues/5).

### Creating a link between arbitrary protein sequences and a genomic reference

If we want to migrate towards a `GRanges` infrastructure, it is
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

## Data
- Availble on prot-main.
- Copy locally in `Pbase/sandbox` and use local path to enable reproducibility.
