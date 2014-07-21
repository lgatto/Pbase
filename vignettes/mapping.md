---
title: "Mapping proteins to their genomic coordinates"
author: "Laurent Gatto - <lg390@cam.ac.uk>"
output:
  html_document:
    toc: true
    theme: united
---



## Introduction

The aim of this document is to document to mapping of proteins and the
tandem mass spectrometry derived peptides they have been inferred from
to genomic locations. 

<img src="figure/schema.png" title="Mapping proteins to a genome reference." alt="Mapping proteins to a genome reference." style="display: block; margin: auto;" />

## Proteins and genome data



We will use a small object from `Pbase` to illustrate how to retrieve
genome coordinates and map a protein back to genomic coordinates. In
the remainder of this vignette, we will concentrate on a spectific
protein, namely P00558.


```r
library("Pbase")
data(p)
p
```

```
## S4 class type     : Proteins
## Class version     : 0.1
## Created           : Wed Jul 16 02:58:33 2014
## Number of Proteins: 9
## Sequences:
##   [1] A4UGR9 [2] A6H8Y1 ... [8] P04075-2 [9] P60709
## Sequence features:
##   [1] DB [2] AccessionNumber ... [10] Comment [11] Filename
## Peptide features:
##   [1] DB [2] AccessionNumber ... [46] spectrumFile [47] databaseFile
```

```r
seqnames(p)
```

```
## [1] "A4UGR9"   "A6H8Y1"   "O43707"   "O75369"   "P00558"   "P02545"  
## [7] "P04075"   "P04075-2" "P60709"
```

```r
dat <- data.frame(i = 1:length(p),
                  npeps = elementLengths(pfeatures(p)),
                  protln = width(aa(p)))
dat
```

```
##          i npeps protln
## A4UGR9   1    36   3374
## A6H8Y1   2    23   2624
## O43707   3     6    911
## O75369   4    13   2602
## P00558   5     5    417
## P02545   6    12    664
## P04075   7    21    364
## P04075-2 8    20    418
## P60709   9     1    375
```

```r
sp <- seqnames(p)[5]
```

### Querying with `biomaRt`

The input information consists of a UniProt identifier and a set of
peptides positions along the protein. The first step is to map the
protein accession number to a gene identifier (here, we will use
Ensembl) and to obtain genome coordinates.

Multiple solutiona are offered to use. Below, we use the `biomaRt`
Bioconductor package.


```r
library("biomaRt")
ens <- useMart("ensembl", "hsapiens_gene_ensembl")

bm <- select(ens, keys = sp,
             keytype = "uniprot_swissprot_accession",
             columns = c(
                 "ensembl_gene_id",
                 "ensembl_transcript_id",
                 "ensembl_peptide_id",
                 "ensembl_exon_id",  
                 "description",     
                 "chromosome_name",
                 "strand",
                 "exon_chrom_start",
                 "exon_chrom_end",
                 "5_utr_start",
                 "5_utr_end",
                 "3_utr_start",
                 "3_utr_end",
                 "start_position",     
                 "end_position"))

bm
```

```
##    ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
## 1  ENSG00000102144       ENST00000373316    ENSP00000362413
## 2  ENSG00000102144       ENST00000373316    ENSP00000362413
## 3  ENSG00000102144       ENST00000373316    ENSP00000362413
## 4  ENSG00000102144       ENST00000373316    ENSP00000362413
## 5  ENSG00000102144       ENST00000373316    ENSP00000362413
## 6  ENSG00000102144       ENST00000373316    ENSP00000362413
## 7  ENSG00000102144       ENST00000373316    ENSP00000362413
## 8  ENSG00000102144       ENST00000373316    ENSP00000362413
## 9  ENSG00000102144       ENST00000373316    ENSP00000362413
## 10 ENSG00000102144       ENST00000373316    ENSP00000362413
## 11 ENSG00000102144       ENST00000373316    ENSP00000362413
## 12 ENSG00000269666       ENST00000597340    ENSP00000472461
## 13 ENSG00000269666       ENST00000597340    ENSP00000472461
## 14 ENSG00000269666       ENST00000597340    ENSP00000472461
## 15 ENSG00000269666       ENST00000597340    ENSP00000472461
## 16 ENSG00000269666       ENST00000597340    ENSP00000472461
## 17 ENSG00000269666       ENST00000597340    ENSP00000472461
## 18 ENSG00000269666       ENST00000597340    ENSP00000472461
## 19 ENSG00000269666       ENST00000597340    ENSP00000472461
## 20 ENSG00000269666       ENST00000597340    ENSP00000472461
## 21 ENSG00000269666       ENST00000597340    ENSP00000472461
## 22 ENSG00000269666       ENST00000597340    ENSP00000472461
##    ensembl_exon_id                                             description
## 1  ENSE00001600900 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 2  ENSE00003506377 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 3  ENSE00003502842 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 4  ENSE00003512377 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 5  ENSE00003581136 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 6  ENSE00003597777 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 7  ENSE00000672997 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 8  ENSE00000672996 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 9  ENSE00000672995 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 10 ENSE00000672994 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 11 ENSE00001948816 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 12 ENSE00003177038 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 13 ENSE00003462979 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 14 ENSE00003487079 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 15 ENSE00003549531 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 16 ENSE00003639123 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 17 ENSE00003515683 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 18 ENSE00003045537 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 19 ENSE00003016554 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 20 ENSE00003003175 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 21 ENSE00003138826 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
## 22 ENSE00003081515 phosphoglycerate kinase 1 [Source:HGNC Symbol;Acc:8896]
##    chromosome_name strand exon_chrom_start exon_chrom_end 5_utr_start
## 1                X      1         77359671       77359902    77359671
## 2                X      1         77365364       77365414          NA
## 3                X      1         77369241       77369396          NA
## 4                X      1         77369513       77369657          NA
## 5                X      1         77372809       77372912          NA
## 6                X      1         77373548       77373667          NA
## 7                X      1         77378332       77378446          NA
## 8                X      1         77378692       77378871          NA
## 9                X      1         77380371       77380548          NA
## 10               X      1         77380824       77380922          NA
## 11               X      1         77381287       77384793          NA
## 12    HG1426_PATCH      1         77365128       77365359    77365128
## 13    HG1426_PATCH      1         77370821       77370871          NA
## 14    HG1426_PATCH      1         77374698       77374853          NA
## 15    HG1426_PATCH      1         77374970       77375114          NA
## 16    HG1426_PATCH      1         77378266       77378369          NA
## 17    HG1426_PATCH      1         77379005       77379124          NA
## 18    HG1426_PATCH      1         77383789       77383903          NA
## 19    HG1426_PATCH      1         77384149       77384328          NA
## 20    HG1426_PATCH      1         77385828       77386005          NA
## 21    HG1426_PATCH      1         77386281       77386379          NA
## 22    HG1426_PATCH      1         77386744       77390250          NA
##    5_utr_end 3_utr_start 3_utr_end start_position end_position
## 1   77359837          NA        NA       77320685     77384793
## 2         NA          NA        NA       77320685     77384793
## 3         NA          NA        NA       77320685     77384793
## 4         NA          NA        NA       77320685     77384793
## 5         NA          NA        NA       77320685     77384793
## 6         NA          NA        NA       77320685     77384793
## 7         NA          NA        NA       77320685     77384793
## 8         NA          NA        NA       77320685     77384793
## 9         NA          NA        NA       77320685     77384793
## 10        NA          NA        NA       77320685     77384793
## 11        NA    77381328  77384793       77320685     77384793
## 12  77365294          NA        NA       77326142     77390250
## 13        NA          NA        NA       77326142     77390250
## 14        NA          NA        NA       77326142     77390250
## 15        NA          NA        NA       77326142     77390250
## 16        NA          NA        NA       77326142     77390250
## 17        NA          NA        NA       77326142     77390250
## 18        NA          NA        NA       77326142     77390250
## 19        NA          NA        NA       77326142     77390250
## 20        NA          NA        NA       77326142     77390250
## 21        NA          NA        NA       77326142     77390250
## 22        NA    77386785  77390250       77326142     77390250
```

We will only consider the gene located on the X chromosome and ignore
the
[sequence](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/info/patches.shtml)
[patch](http://www.ncbi.nlm.nih.gov/nuccore/NW_003871101.3?report=genbank)
and conveniently store the transcript and protein identifiers for
later use.


```r
bm <- bm[bm$chromosome_name == "X", ]
tr <- unique(bm$ensembl_transcript_id)
pr <- unique(bm$ensembl_peptide_id)
```

The information retrieved gives us various information, in particular
the Ensembl gene identifer ENSG00000102144 and its starting
and ending position 77320685 and 77384793
on chromosome X.

### Genome visualisation with `Gviz`

While the above defines all necessary genomic coordinates, and as we
will make use of the `Gviz` plotting infrastructure, it is more
straightforward to fetch the above information via `biomaRt` using
specialised `Gviz` infrastructure.


```r
library("Gviz")
options(ucscChromosomeNames=FALSE)

bmTrack <- BiomartGeneRegionTrack(start = min(bm$start_position),
                                  end = max(bm$end_position),
                                  chromosome = "chrX",
                                  genome = "hg19")

bmTracks <- split(bmTrack, transcript(bmTrack))

grTrack <- bmTracks[[which(names(bmTracks) == tr)]]
names(grTrack) <- tr
```

We can convince ourselves that the manual retrieved data in `bm` and
the data in the `BiomartRegionTrack` subsequently subset into a
`GeneRegionTrack` match.


```r
as(grTrack, "data.frame")
```

```
##    X.seqnames  X.start    X.end X.width X.strand      X.feature
## 1           X 77359671 77359837     167        +           utr5
## 2           X 77359838 77359902      65        + protein_coding
## 3           X 77365364 77365414      51        + protein_coding
## 4           X 77369241 77369396     156        + protein_coding
## 5           X 77369513 77369657     145        + protein_coding
## 6           X 77372809 77372912     104        + protein_coding
## 7           X 77373548 77373667     120        + protein_coding
## 8           X 77378332 77378446     115        + protein_coding
## 9           X 77378692 77378871     180        + protein_coding
## 10          X 77380371 77380548     178        + protein_coding
## 11          X 77380824 77380922      99        + protein_coding
## 12          X 77381287 77381327      41        + protein_coding
## 13          X 77381328 77384793    3466        +           utr3
##             X.gene          X.exon    X.transcript X.symbol X.rank X.phase
## 1  ENSG00000102144 ENSE00001600900 ENST00000373316     PGK1      1      -1
## 2  ENSG00000102144 ENSE00001600900 ENST00000373316     PGK1      1      -1
## 3  ENSG00000102144 ENSE00003506377 ENST00000373316     PGK1      2       2
## 4  ENSG00000102144 ENSE00003502842 ENST00000373316     PGK1      3       2
## 5  ENSG00000102144 ENSE00003512377 ENST00000373316     PGK1      4       2
## 6  ENSG00000102144 ENSE00003581136 ENST00000373316     PGK1      5       0
## 7  ENSG00000102144 ENSE00003597777 ENST00000373316     PGK1      6       2
## 8  ENSG00000102144 ENSE00000672997 ENST00000373316     PGK1      7       2
## 9  ENSG00000102144 ENSE00000672996 ENST00000373316     PGK1      8       0
## 10 ENSG00000102144 ENSE00000672995 ENST00000373316     PGK1      9       0
## 11 ENSG00000102144 ENSE00000672994 ENST00000373316     PGK1     10       1
## 12 ENSG00000102144 ENSE00001948816 ENST00000373316     PGK1     11       0
## 13 ENSG00000102144 ENSE00001948816 ENST00000373316     PGK1     11      -1
```

Below, we create additional ideogram and axis tracks and produce a
visualisation of our genomic coordinates of interest.


```r
ideoTrack <- IdeogramTrack(genome = "hg19",
                           chromosome = "chrX")
axisTrack <- GenomeAxisTrack()
seqlevels(ranges(grTrack)) <- 
    chromosome(grTrack) <-
        chromosome(ideoTrack)
plotTracks(list(ideoTrack, axisTrack, grTrack),
           add53 = TRUE, add35 = TRUE)
```

<img src="figure/gvizfig.png" title="plot of chunk gvizfig" alt="plot of chunk gvizfig" style="display: block; margin: auto;" />

### Using `TranscriptDb` instances

Finally, on could also use `TranscriptDb` objects to create the
`GeneRegionTrack`. Below, we use the UCSC annotation (which is
directly available from Bioconductor; Ensembl TranscriptDb instances
could be generated manually - **TODO** document this) to extract our
regions of interest.


```r
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txTr <- GeneRegionTrack(txdb, chromosome = "chrX",
                        start = bm$start_position[1],
                        end = bm$end_position[1])

head(as(txTr, "data.frame"))
```

```
##   X.seqnames  X.start    X.end X.width X.strand X.feature    X.id
## 1       chrX 77359666 77359837     172        +      utr5 unknown
## 2       chrX 77359838 77359902      65        +       CDS unknown
## 3       chrX 77361859 77362168     310        +      utr5 unknown
## 4       chrX 77365364 77365382      19        +      utr5 unknown
## 5       chrX 77365364 77365414      51        +       CDS unknown
## 6       chrX 77365383 77365414      32        +       CDS unknown
##         X.exon X.transcript X.gene   X.symbol X.density
## 1 uc004ecz.4_1   uc004ecz.4   5230 uc004ecz.4         1
## 2 uc004ecz.4_1   uc004ecz.4   5230 uc004ecz.4         1
## 3 uc011mqq.2_1   uc011mqq.2   5230 uc011mqq.2         1
## 4 uc011mqq.2_2   uc011mqq.2   5230 uc011mqq.2         1
## 5 uc004ecz.4_2   uc004ecz.4   5230 uc004ecz.4         1
## 6 uc011mqq.2_2   uc011mqq.2   5230 uc011mqq.2         1
```

## Peptide features



The peptides that have been experimentally observed are available as
ranges (coordinates) along the protein sequences. For example, below,
we see 5 peptides (ELNYFAKALESPER, DLMSKAEK, QIVWNGPVGVFEWEAFAR, FHVEEEGKGKDASGNK, GTKALMDEVVK) have been identified for our
protein of interest P00558.


```r
pp <- p[sp]
pranges(pp)
```

```
## IRangesList of length 1
## $P00558
## IRanges of length 5
##     start end width  names
## [1]   193 206    14 P00558
## [2]   268 275     8 P00558
## [3]   333 350    18 P00558
## [4]   124 139    16 P00558
## [5]   351 361    11 P00558
```

```r
plot(pp)
```

<img src="figure/pepsfig.png" title="plot of chunk pepsfig" alt="plot of chunk pepsfig" style="display: block; margin: auto;" />

The aim of this document is to document the mapping of peptides,
i.e. ranges along a protein sequence to ranges along the genome
reference. In other words, our aim is the convert protein coordinates
to genome coordinates.

## Mapping peptides back to the genome

Additional features of interest (in our case peptides) can be added to
the genomic visualisation using dedicated annotation tracks, that can
be constructed as shown below. 

### Comparing protein and translated DNA sequences

The first check that we want to implement is to verify that we can
regenerate the protein amino acid sequence from the genome regions
that we have extracted. We start by subsetting only the actual protein
coding regions, i.e ignoring the 5' and 3' untranslated regions of our
Genome Region track.


```r
(prng <- ranges(grTrack[feature(grTrack) == "protein_coding",]))
```

```
## GRanges with 11 ranges and 7 metadata columns:
##        seqnames               ranges strand   |        feature
##           <Rle>            <IRanges>  <Rle>   |    <character>
##    [1]     chrX [77359838, 77359902]      +   | protein_coding
##    [2]     chrX [77365364, 77365414]      +   | protein_coding
##    [3]     chrX [77369241, 77369396]      +   | protein_coding
##    [4]     chrX [77369513, 77369657]      +   | protein_coding
##    [5]     chrX [77372809, 77372912]      +   | protein_coding
##    ...      ...                  ...    ... ...            ...
##    [7]     chrX [77378332, 77378446]      +   | protein_coding
##    [8]     chrX [77378692, 77378871]      +   | protein_coding
##    [9]     chrX [77380371, 77380548]      +   | protein_coding
##   [10]     chrX [77380824, 77380922]      +   | protein_coding
##   [11]     chrX [77381287, 77381327]      +   | protein_coding
##                   gene            exon      transcript      symbol
##            <character>     <character>     <character> <character>
##    [1] ENSG00000102144 ENSE00001600900 ENST00000373316        PGK1
##    [2] ENSG00000102144 ENSE00003506377 ENST00000373316        PGK1
##    [3] ENSG00000102144 ENSE00003502842 ENST00000373316        PGK1
##    [4] ENSG00000102144 ENSE00003512377 ENST00000373316        PGK1
##    [5] ENSG00000102144 ENSE00003581136 ENST00000373316        PGK1
##    ...             ...             ...             ...         ...
##    [7] ENSG00000102144 ENSE00000672997 ENST00000373316        PGK1
##    [8] ENSG00000102144 ENSE00000672996 ENST00000373316        PGK1
##    [9] ENSG00000102144 ENSE00000672995 ENST00000373316        PGK1
##   [10] ENSG00000102144 ENSE00000672994 ENST00000373316        PGK1
##   [11] ENSG00000102144 ENSE00001948816 ENST00000373316        PGK1
##             rank     phase
##        <numeric> <integer>
##    [1]         1        -1
##    [2]         2         2
##    [3]         3         2
##    [4]         4         2
##    [5]         5         0
##    ...       ...       ...
##    [7]         7         2
##    [8]         8         0
##    [9]         9         0
##   [10]        10         1
##   [11]        11         0
##   ---
##   seqlengths:
##    chrX
##      NA
```

We also need the actual genome sequence (so far, we have only dealt
with regions and features). There are multiple ways to obtain genome
sequences with R and Bioconductor. As we have been working with the
Ensembl reference genome, we could simply
[download](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/)
the whole genome of only the chromosome X and load it as an
`AAStringSet` and extract the protein coding regions of interest at
the coordinates defined by the ranges in `prng`.


```r
library("Biostrings")
chrx <- readDNAStringSet("./Homo_sapiens.GRCh37.75.dna.chromosome.X.fa.gz")
seq <- extractAt(chrx[[1]], ranges(prng))
pptr <- translate(unlist(seq))
```

It is also possible to use readily package genomes to do this. Above,
we have use the `TxDb.Hsapiens.UCSC.hg19.knownGene` package defining
transcripts. There is a similar package for the human UCSC genome
build, namely `BSgenome.Hsapiens.UCSC.hg19`. However, as we have
focused on Ensembl data annotation, we will first need to convert our
Ensembl transcript (or protein) identifier ENST00000373316 (ENSP00000362413) into the
UCSC equivalent. (**TODO** Use also `Homo.sapiens` annotation package.)


```r
ucsc <- select(ens, key = pr,
               keytype = "ensembl_peptide_id",
               columns = "ucsc")
ucsc
```

```
##         ucsc
## 1 uc004ecz.4
```

```r
library("BSgenome.Hsapiens.UCSC.hg19")
prng2 <- ranges(txTr[transcript(txTr) %in% ucsc])
prng2 <- prng2[prng2$feature == "CDS"]

seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, prng2)
pptr <- translate(unlist(seq))
## remove stop codon
pptr <- pptr[1:417]
```


```r
writePairwiseAlignments(pairwiseAlignment(pp[[1]], pptr))
```

```
## ########################################
## # Program: Biostrings (version 2.33.12), a Bioconductor package
## # Rundate: Mon Jul 21 01:43:05 2014
## ########################################
## #=======================================
## #
## # Aligned_sequences: 2
## # 1: P1
## # 2: S1
## # Matrix: NA
## # Gap_penalty: 14.0
## # Extend_penalty: 4.0
## #
## # Length: 417
## # Identity:     417/417 (100.0%)
## # Similarity:    NA/417 (NA%)
## # Gaps:           0/417 (0.0%)
## # Score: 1795
## #
## #
## #=======================================
## 
## P1                 1 MSLSNKLTLDKLDVKGKRVVMRVDFNVPMKNNQITNNQRIKAAVPSIKFC     50
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1                 1 MSLSNKLTLDKLDVKGKRVVMRVDFNVPMKNNQITNNQRIKAAVPSIKFC     50
## 
## P1                51 LDNGAKSVVLMSHLGRPDGVPMPDKYSLEPVAVELKSLLGKDVLFLKDCV    100
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1                51 LDNGAKSVVLMSHLGRPDGVPMPDKYSLEPVAVELKSLLGKDVLFLKDCV    100
## 
## P1               101 GPEVEKACANPAAGSVILLENLRFHVEEEGKGKDASGNKVKAEPAKIEAF    150
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1               101 GPEVEKACANPAAGSVILLENLRFHVEEEGKGKDASGNKVKAEPAKIEAF    150
## 
## P1               151 RASLSKLGDVYVNDAFGTAHRAHSSMVGVNLPQKAGGFLMKKELNYFAKA    200
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1               151 RASLSKLGDVYVNDAFGTAHRAHSSMVGVNLPQKAGGFLMKKELNYFAKA    200
## 
## P1               201 LESPERPFLAILGGAKVADKIQLINNMLDKVNEMIIGGGMAFTFLKVLNN    250
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1               201 LESPERPFLAILGGAKVADKIQLINNMLDKVNEMIIGGGMAFTFLKVLNN    250
## 
## P1               251 MEIGTSLFDEEGAKIVKDLMSKAEKNGVKITLPVDFVTADKFDENAKTGQ    300
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1               251 MEIGTSLFDEEGAKIVKDLMSKAEKNGVKITLPVDFVTADKFDENAKTGQ    300
## 
## P1               301 ATVASGIPAGWMGLDCGPESSKKYAEAVTRAKQIVWNGPVGVFEWEAFAR    350
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1               301 ATVASGIPAGWMGLDCGPESSKKYAEAVTRAKQIVWNGPVGVFEWEAFAR    350
## 
## P1               351 GTKALMDEVVKATSRGCITIIGGGDTATCCAKWNTEDKVSHVSTGGGASL    400
##                      ||||||||||||||||||||||||||||||||||||||||||||||||||
## S1               351 GTKALMDEVVKATSRGCITIIGGGDTATCCAKWNTEDKVSHVSTGGGASL    400
## 
## P1               401 ELLEGKVLPGVDALSNI    417
##                      |||||||||||||||||
## S1               401 ELLEGKVLPGVDALSNI    417
## 
## 
## #---------------------------------------
## #---------------------------------------
```

### Calculating new coordinates

To calculate the genomic coordinates of peptides we

1) Calculate the contiguous exon coordinates on the cDNA sequence
   `prex`.

2) Calculate the peptide corrdinates on the cDNA sequence `peprgn2`
   from the original peptide ranges along the protein sequence.

3) We calculate the indices of the exons in which all the peptides
   start and end on the cDNA sequence `start_ex` and `end_ex`.

4) We infer the coordinates on the genome from the actual peptide
   starting/ending position and exon index on the cDNA sequence, the
   exons coordinates on the cDNA sequences and the matching exons
   coordinates on the genome and store these in an `IRanges` object.


```r
## exons ranges along the protein
j <- cumsum(width(seq))
i <- cumsum(c(1, width(seq)))[1:length(j)]
prex <- IRanges(start = i, end = j)

peprng <- pranges(pp)[[1]]
## peptide positions on cdna
peprng2 <- IRanges(start = 1 + (start(peprng)-1) * 3,
                   width = width(peprng) * 3)

## find exon and position in prex 
start_ex <- subjectHits(findOverlaps(start(peprng2), prex))
end_ex <- subjectHits(findOverlaps(end(peprng2), prex))

getPos <- function(p, i, prtex = prex, nclex = prng2) {
    ## position in cdna
    ## exon index
    start(prng2[i]) + (p - start(prtex[i]))
}

peptides_on_genome <- IRanges(start = getPos(start(peprng2), start_ex),
                              end = getPos(end(peprng2), end_ex))
```

### Plotting

Based on the new peptide genomic coordinates, it is now
straightforward to create a new `AnnotationTrack` and add it the the
track visualisation.


```r
pepTr <- AnnotationTrack(start = start(peptides_on_genome),
                         end = end(peptides_on_genome),
                         chr = "chrX", genome = "hg19",
                         strand = "*",
                         id = pcols(pp)[[1]]$pepseq,
                         name = "pfeatures",
                         col = "steelblue")


plotTracks(list(ideoTrack, axisTrack, grTrack, pepTr),
           groupAnnotation = "id",
           just.group = "below",
           fontsize.group = 9,
           add53 = TRUE, add35 = TRUE)
```

<img src="figure/pepcoords.png" title="plot of chunk pepcoords" alt="plot of chunk pepcoords" style="display: block; margin: auto;" />

### DetailsAnnotationTrack

Finally, we customise the figure by adding a track with the $MS^2$
spectra. The raw data used to search the protein database an create
`p` are available as an `MSnExp` object. 


```r
data(pms)

library("ggplot2")
details <- function(identifier, ...) {
    p <- plot(pms[[as.numeric(identifier)]], full=TRUE, plot=FALSE) + ggtitle("") 
    p <- p + theme_bw() + theme(axis.text.y = element_blank(),
                                axis.text.x = element_blank()) + 
                                labs(x = NULL, y = NULL)
                                        
    print(p, newpage=FALSE)
}

deTrack <- AnnotationTrack(start = start(peptides_on_genome),
                           end = end(peptides_on_genome),
                           genome = "hg19", chromosom = "chrX",
                           id = pcols(pp)[[1]]$acquisitionnum,
                           name = "MS2 spectra",
                           stacking = "squish", fun = details)

plotTracks(list(ideoTrack, axisTrack, deTrack, grTrack),
           add53 = TRUE, add35 = TRUE)
```

<img src="figure/msmsspectra.png" title="plot of chunk msmsspectra" alt="plot of chunk msmsspectra" style="display: block; margin: auto;" />

**TODO** Check spectra. Describe how data tracks can be used to
  overlay additional information, such as quantitation data,
  identification scores, coverage, ...

## Session information


```r
sessionInfo()
```

```
## R Under development (unstable) (2014-04-10 r65396)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  grid      methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] ggplot2_1.0.0                     BSgenome.Hsapiens.UCSC.hg19_1.4.0
##  [3] BSgenome_1.33.8                   rtracklayer_1.25.13              
##  [5] biomaRt_2.21.1                    Biostrings_2.33.12               
##  [7] XVector_0.5.7                     Pbase_0.1.6                      
##  [9] Rcpp_0.11.2                       Gviz_1.9.10                      
## [11] GenomicRanges_1.17.23             GenomeInfoDb_1.1.12              
## [13] IRanges_1.99.22                   S4Vectors_0.1.2                  
## [15] BiocGenerics_0.11.3               rmarkdown_0.2.46                 
## [17] knitr_1.6                        
## 
## loaded via a namespace (and not attached):
##  [1] affy_1.43.3               affyio_1.33.0            
##  [3] AnnotationDbi_1.27.8      BatchJobs_1.3            
##  [5] BBmisc_1.7                Biobase_2.25.0           
##  [7] BiocInstaller_1.15.5      BiocParallel_0.7.7       
##  [9] biovizBase_1.13.8         bitops_1.0-6             
## [11] brew_1.0-6                checkmate_1.1            
## [13] cleaver_1.3.7             cluster_1.15.2           
## [15] codetools_0.2-8           colorspace_1.2-4         
## [17] data.table_1.9.2          DBI_0.2-7                
## [19] dichromat_2.0-0           digest_0.6.4             
## [21] doParallel_1.0.8          evaluate_0.5.5           
## [23] fail_1.2                  foreach_1.4.2            
## [25] formatR_0.10              Formula_1.1-1            
## [27] GenomicAlignments_1.1.20  GenomicFeatures_1.17.12  
## [29] gtable_0.1.2              Hmisc_3.14-4             
## [31] htmltools_0.2.4           impute_1.39.0            
## [33] iterators_1.0.7           labeling_0.2             
## [35] lattice_0.20-29           latticeExtra_0.6-26      
## [37] limma_3.21.10             MALDIquant_1.10          
## [39] MASS_7.3-33               matrixStats_0.10.0       
## [41] MSnbase_1.13.12           munsell_0.4.2            
## [43] mzID_1.3.2                mzR_1.11.10              
## [45] pcaMethods_1.55.0         plyr_1.8.1               
## [47] preprocessCore_1.27.1     proto_0.3-10             
## [49] Pviz_0.99.0               R.methodsS3_1.6.1        
## [51] RColorBrewer_1.0-5        RCurl_1.95-4.1           
## [53] reshape2_1.4              Rsamtools_1.17.31        
## [55] RSQLite_0.11.4            scales_0.2.4             
## [57] sendmailR_1.1-2           splines_3.2.0            
## [59] stats4_3.2.0              stringr_0.6.2            
## [61] survival_2.37-7           tools_3.2.0              
## [63] VariantAnnotation_1.11.16 vsn_3.33.0               
## [65] XML_3.98-1.1              zlibbioc_1.11.1
```
