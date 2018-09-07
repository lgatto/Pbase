Pbase
=====

Manipulating and exploring protein and proteomics data. 

## Installation

It is advised to install `Pbase` from Bioconductor:

    library("BiocManager")
	BiocManager::install("Pbase")

Note however that you will need the *devel* version of Bioconductor
for this. See `?BiocManager::install` for details.

From github using `devtools::install_github`:

    library("devtools")
    install_github("ComputationalProteomicsUnit/Pbase")

### Dependencies

See the `DESCRIPTION` file for a complete list.

## Getting started

Currently, the best way to get started is `?Proteins` and the
[`Pbase-data`](http://bioconductor.org/packages/devel/bioc/vignettes/Pbase/inst/doc/Pbase-data.html)
vignette. More documentation is on its way.

## Development

`Pbase` is under heavy development and is likely to considerably
change in the near future. Suggestion and bug reports are welcome and
can be filed as
[github issues](https://github.com/ComputationalProteomicsUnit/pbase/issues).

If you would like to contribute, please directly send pull requests
for minor contributions and typos. For major contributions, we suggest
to first get in touch with the package maintainers. 

## Ideas

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

### Protein domains

Available through the integration with the `EnsmbleDb` package. See the `Pbase-with-ensembldb` vignette.

### Mapping a Protein Sequence to a Genome Sequence

See the [`mapping`](http://bioconductor.org/packages/devel/bioc/vignettes/Pbase/inst/doc/mapping.html) vignette.

See also
[this document](https://github.com/ComputationalProteomicsUnit/Intro-Integ-Omics-Prot/blob/master/mapping.md)
for additional examples and integration with RNA-seq data.

## Interoperability

The package allows to easily interact with `AAString` and
`AAStringSet` instances, protein databases such as UniProt (and
possibly biomaRt in the future) using protein identifiers, protein
identification results (`mzID` or (devel) `mzR` packages) and possibly
also `MSnExp` and `MSnSet` instances.

