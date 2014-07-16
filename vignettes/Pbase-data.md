---
title: "Pbase example data"
author: "Laurent Gatto - <lg390@cam.ac.uk>"
output:
  html_document:
    toc: true
    theme: united
---



The original data is a 10 fmol
[Peptide Retention Time Calibration Mixture](http://www.piercenet.com/product/peptide-retention-time-calibration-mixture)
spiked into 50 ng HeLa background acquired on a Thermo Orbitrap Q
Exactive instrument. A restricted set of high scoring human proteins
from the UniProt release 2014_06 were searched using the MSGF+ engine.

## The fasta database


```r
library("Biostrings")
```

```
## Loading required package: XVector
```

```r
fafile <- system.file("extdata/HUMAN_2014_06-09prots.fasta", package = "Pbase")
fa <- readAAStringSet(fafile)
fa
```

```
##   A AAStringSet instance of length 9
##     width seq                                          names               
## [1]   375 MDDDIAALVVDNGSGMCKAGF...WISKQEYDESGPSIVHRKCF sp|P60709|ACTB_HU...
## [2]   664 METPSQRRATRSGAQASSTPL...SYLLGNSSPRTQSPQNCSIM sp|P02545|LMNA_HU...
## [3]   417 MSLSNKLTLDKLDVKGKRVVM...ASLELLEGKVLPGVDALSNI sp|P00558|PGK1_HU...
## [4]  2602 MPVTEKDLAEDAPWKKIQQNT...LAVKWGEEHIPGSPFHVTVP sp|O75369|FLNB_HU...
## [5]   911 MVDYHAANQSYQYGPSSAGNG...VPGALDYKSFSTALYGESDL sp|O43707|ACTN4_H...
## [6]  2624 MFRRARLSVKPNVRPGVGARG...ATTVSEYFFNDIFIEVDETE sp|A6H8Y1|BDP1_HU...
## [7]  3374 MSPESGHSRIFEATAGPNKPE...TLSKDSLSNGVPSGRQAEFS sp|A4UGR9|XIRP2_H...
## [8]   364 MPYQYPALTPEQKKELSDIAH...PSGQAGAAASESLFVSNHAY sp|P04075|ALDOA_H...
## [9]   418 MARRKPEGSSFNMTHLSMAMA...PSGQAGAAASESLFVSNHAY sp|P04075-2|ALDOA...
```

## The PSM data


```r
library("mzID")
idfile <- system.file("extdata/Thermo_Hela_PRTC_1_selected.mzid", package = "Pbase")
id <- flatten(mzID(idfile))
```

```
## reading Thermo_Hela_PRTC_1_selected.mzid... DONE!
```

```r
dim(id)
```

```
## [1] 137  39
```

```r
head(id)
```

```
##   spectrumid acquisitionnum passthreshold rank calculatedmasstocharge
## 1  index=173            173          TRUE    1      1136.574462890625
## 2  index=173            173          TRUE    1      1136.574462890625
## 3  index=163            163          TRUE    1      1136.574462890625
## 4  index=163            163          TRUE    1      1136.574462890625
## 5    index=8              8          TRUE    1      469.2348937988281
## 6    index=8              8          TRUE    1      469.2348937988281
##   experimentalmasstocharge chargestate ms-gf:denovoscore ms-gf:evalue
## 1        1137.066650390625           2               132    2.597e-18
## 2        1137.066650390625           2               132    2.597e-18
## 3        1136.574462890625           2               230    4.943e-17
## 4        1136.574462890625           2               230    4.943e-17
## 5       469.23480224609375           2                49    3.824e-06
## 6       469.23480224609375           2                49    3.824e-06
##   ms-gf:rawscore ms-gf:specevalue assumeddissociationmethod
## 1            118        2.277e-22                       CID
## 2            118        2.277e-22                       CID
## 3            186        4.333e-21                       CID
## 4            186        4.333e-21                       CID
## 5             44        3.374e-10                       CID
## 6             44        3.374e-10                       CID
##   ctermioncurrentratio explainedioncurrentratio isotopeerror meanerrorall
## 1            0.5061673                0.5888804            1    17.757027
## 2            0.5061673                0.5888804            1    17.757027
## 3           0.48202217               0.57532895            0    18.029726
## 4           0.48202217               0.57532895            0    18.029726
## 5             0.285157               0.44413796            0    4.1399555
## 6             0.285157               0.44413796            0    4.1399555
##   meanerrortop7 meanrelerrorall meanrelerrortop7 ms2ioncurrent
## 1     4.3787823     -15.7646885        -3.255431      526864.0
## 2     4.3787823     -15.7646885        -3.255431      526864.0
## 3     5.0083427       -16.59046        -3.333969     3342169.2
## 4     5.0083427       -16.59046        -3.333969     3342169.2
## 5      4.087293      -2.5620477       -2.0585544     8156462.0
## 6      4.087293      -2.5620477       -2.0585544     8156462.0
##   ntermioncurrentratio nummatchedmainions stdeverrorall stdeverrortop7
## 1           0.08271315                 25     23.135044       3.582861
## 2           0.08271315                 25     23.135044       3.582861
## 3           0.09330681                 30     28.647552      2.9929223
## 4           0.09330681                 30     28.647552      2.9929223
## 5           0.15898097                  9     2.7658114      2.7357001
## 6           0.15898097                  9     2.7658114      2.7357001
##   stdevrelerrorall stdevrelerrortop7                  pepseq modified
## 1         24.53603          4.627396 GVVPLAGTNGETTTQGLDGLSER    FALSE
## 2         24.53603          4.627396 GVVPLAGTNGETTTQGLDGLSER    FALSE
## 3        29.504406          4.788082 GVVPLAGTNGETTTQGLDGLSER    FALSE
## 4        29.504406          4.788082 GVVPLAGTNGETTTQGLDGLSER    FALSE
## 5         4.269058          4.466808                AAQEEYVK    FALSE
## 6         4.269058          4.466808                AAQEEYVK    FALSE
##   modification isdecoy post pre end start               accession length
## 1         <NA>   FALSE    C   K 134   112   sp|P04075|ALDOA_HUMAN    364
## 2         <NA>   FALSE    C   K 188   166 sp|P04075-2|ALDOA_HUMAN    418
## 3         <NA>   FALSE    C   K 134   112   sp|P04075|ALDOA_HUMAN    364
## 4         <NA>   FALSE    C   K 188   166 sp|P04075-2|ALDOA_HUMAN    418
## 5         <NA>   FALSE    R   K 330   323   sp|P04075|ALDOA_HUMAN    364
## 6         <NA>   FALSE    R   K 384   377 sp|P04075-2|ALDOA_HUMAN    418
##                                                              description
## 1    Fructose-bisphosphate aldolase A OS=Homo sapiens GN=ALDOA PE=1 SV=2
## 2 Isoform 2 of Fructose-bisphosphate aldolase A OS=Homo sapiens GN=ALDOA
## 3    Fructose-bisphosphate aldolase A OS=Homo sapiens GN=ALDOA PE=1 SV=2
## 4 Isoform 2 of Fructose-bisphosphate aldolase A OS=Homo sapiens GN=ALDOA
## 5    Fructose-bisphosphate aldolase A OS=Homo sapiens GN=ALDOA PE=1 SV=2
## 6 Isoform 2 of Fructose-bisphosphate aldolase A OS=Homo sapiens GN=ALDOA
##                    spectrumFile                databaseFile
## 1 Thermo_Hela_PRTC_selected.mgf HUMAN_2014_06-09prots.fasta
## 2 Thermo_Hela_PRTC_selected.mgf HUMAN_2014_06-09prots.fasta
## 3 Thermo_Hela_PRTC_selected.mgf HUMAN_2014_06-09prots.fasta
## 4 Thermo_Hela_PRTC_selected.mgf HUMAN_2014_06-09prots.fasta
## 5 Thermo_Hela_PRTC_selected.mgf HUMAN_2014_06-09prots.fasta
## 6 Thermo_Hela_PRTC_selected.mgf HUMAN_2014_06-09prots.fasta
```

## The Proteins object


```r
library("Pbase")
p <- Proteins(fafile)
p <- addIdentificationData(p, idfile)
```

```
## reading Thermo_Hela_PRTC_1_selected.mzid... DONE!
```

```r
p
```

```
## S4 class type     : Proteins
## Class version     : 0.1
## Created           : Wed Jul 16 02:45:51 2014
## Number of Proteins: 9
## Sequences:
##   [1] A4UGR9 [2] A6H8Y1 ... [8] P04075-2 [9] P60709
## Sequence features:
##   [1] DB [2] AccessionNumber ... [10] Comment [11] Filename
## Peptide features:
##   [1] DB [2] AccessionNumber ... [46] spectrumFile [47] databaseFile
```

A `Proteins` object is composed of a set of protein sequences
accessible with the `aa` accessor as well as an optional set of
peptides features that are mapped as coordinates along the proteins,
available with `pranges`. The actual peptide sequences can be extraced
with `pfeatures`.


```r
aa(p)
```

```
##   A AAStringSet instance of length 9
##     width seq                                          names               
## [1]  3374 MSPESGHSRIFEATAGPNKPE...TLSKDSLSNGVPSGRQAEFS A4UGR9
## [2]  2624 MFRRARLSVKPNVRPGVGARG...ATTVSEYFFNDIFIEVDETE A6H8Y1
## [3]   911 MVDYHAANQSYQYGPSSAGNG...VPGALDYKSFSTALYGESDL O43707
## [4]  2602 MPVTEKDLAEDAPWKKIQQNT...LAVKWGEEHIPGSPFHVTVP O75369
## [5]   417 MSLSNKLTLDKLDVKGKRVVM...ASLELLEGKVLPGVDALSNI P00558
## [6]   664 METPSQRRATRSGAQASSTPL...SYLLGNSSPRTQSPQNCSIM P02545
## [7]   364 MPYQYPALTPEQKKELSDIAH...PSGQAGAAASESLFVSNHAY P04075
## [8]   418 MARRKPEGSSFNMTHLSMAMA...PSGQAGAAASESLFVSNHAY P04075-2
## [9]   375 MDDDIAALVVDNGSGMCKAGF...WISKQEYDESGPSIVHRKCF P60709
```

```r
pranges(p)
```

```
## IRangesList of length 9
## $A4UGR9
## IRanges of length 36
##      start  end width  names
## [1]   2743 2760    18 A4UGR9
## [2]    307  318    12 A4UGR9
## [3]   1858 1870    13 A4UGR9
## [4]   1858 1870    13 A4UGR9
## [5]   1699 1708    10 A4UGR9
## ...    ...  ...   ...    ...
## [32]  2082 2094    13 A4UGR9
## [33]    20   31    12 A4UGR9
## [34]  1712 1729    18 A4UGR9
## [35]    48   61    14 A4UGR9
## [36]  2743 2756    14 A4UGR9
## 
## ...
## <8 more elements>
```

```r
pfeatures(p)
```

```
## AAStringSetList of length 9
## [["A4UGR9"]] QEITQNKSFFSSVKESQR LPVPKDVYSKQR ... QEITQNKSFFSSVK
## [["A6H8Y1"]] EDAEQVALEVDLNQKKRR ... ARLSVKPNVRPGVGARGSTASNPQRGR
## [["O43707"]] QQRKTFTAWCNSHLR CQKICDQWDALGSLTHSR ... VGWEQLLTTIAR
## [["O75369"]] DLDIIDNYDYSHTVK PFDLVIPFAVRK ... VQAQGPGLKEAFTNK
## [["P00558"]] ELNYFAKALESPER DLMSKAEK QIVWNGPVGVFEWEAFAR FHVEEEGKGKDASGNK GTKALMDEVVK
## [["P02545"]] METPSQRRATR DTSRRLLAEKEREMAEMR ... RATRSGAQASSTPLSPTR
## [["P04075"]] GVVPLAGTNGETTTQGLDGLSER ... TVPPAVTGITFLSGGQSEEEASINLNAINK
## [["P04075-2"]] GVVPLAGTNGETTTQGLDGLSER ... TVPPAVTGITFLSGGQSEEEASINLNAINK
## [["P60709"]] DLTDYLMKILTER
```

A Proteins instance is further described by general
`metadata`. Protein sequence and peptide features annotations can be
accessed with `ametadata` and `pmetadata` (or `acols` and `pcols`)
respectively.


```r
metadata(p)
```

```
## $created
## [1] "Wed Jul 16 02:45:51 2014"
```

```r
head(acols(p))
```

```
## DataFrame with 6 rows and 11 columns
##      DB AccessionNumber   EntryName IsoformName
##   <Rle>     <character> <character>       <Rle>
## 1    sp          A4UGR9 XIRP2_HUMAN          NA
## 2    sp          A6H8Y1  BDP1_HUMAN          NA
## 3    sp          O43707 ACTN4_HUMAN          NA
## 4    sp          O75369  FLNB_HUMAN          NA
## 5    sp          P00558  PGK1_HUMAN          NA
## 6    sp          P02545  LMNA_HUMAN          NA
##                                         ProteinName OrganismName GeneName
##                                         <character>        <Rle>    <Rle>
## 1     Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 2 Transcription factor TFIIIB component B'' homolog Homo sapiens     BDP1
## 3                                   Alpha-actinin-4 Homo sapiens    ACTN4
## 4                                         Filamin-B Homo sapiens     FLNB
## 5                         Phosphoglycerate kinase 1 Homo sapiens     PGK1
## 6                                      Prelamin-A/C Homo sapiens     LMNA
##            ProteinExistence SequenceVersion Comment
##                       <Rle>           <Rle>   <Rle>
## 1 Evidence at protein level               2      NA
## 2 Evidence at protein level               3      NA
## 3 Evidence at protein level               2      NA
## 4 Evidence at protein level               2      NA
## 5 Evidence at protein level               3      NA
## 6 Evidence at protein level               1      NA
##                                                                                        Filename
##                                                                                           <Rle>
## 1 /home/lgatto/R/x86_64-unknown-linux-gnu-library/3.2/Pbase/extdata/HUMAN_2014_06-09prots.fasta
## 2 /home/lgatto/R/x86_64-unknown-linux-gnu-library/3.2/Pbase/extdata/HUMAN_2014_06-09prots.fasta
## 3 /home/lgatto/R/x86_64-unknown-linux-gnu-library/3.2/Pbase/extdata/HUMAN_2014_06-09prots.fasta
## 4 /home/lgatto/R/x86_64-unknown-linux-gnu-library/3.2/Pbase/extdata/HUMAN_2014_06-09prots.fasta
## 5 /home/lgatto/R/x86_64-unknown-linux-gnu-library/3.2/Pbase/extdata/HUMAN_2014_06-09prots.fasta
## 6 /home/lgatto/R/x86_64-unknown-linux-gnu-library/3.2/Pbase/extdata/HUMAN_2014_06-09prots.fasta
```

```r
head(pcols(p))
```

```
## SplitDataFrameList of length 6
## $A4UGR9
## DataFrame with 36 rows and 47 columns
##        DB AccessionNumber   EntryName IsoformName
##     <Rle>     <character> <character>       <Rle>
## 1      sp          A4UGR9 XIRP2_HUMAN          NA
## 2      sp          A4UGR9 XIRP2_HUMAN          NA
## 3      sp          A4UGR9 XIRP2_HUMAN          NA
## 4      sp          A4UGR9 XIRP2_HUMAN          NA
## 5      sp          A4UGR9 XIRP2_HUMAN          NA
## ...   ...             ...         ...         ...
## 32     sp          A4UGR9 XIRP2_HUMAN          NA
## 33     sp          A4UGR9 XIRP2_HUMAN          NA
## 34     sp          A4UGR9 XIRP2_HUMAN          NA
## 35     sp          A4UGR9 XIRP2_HUMAN          NA
## 36     sp          A4UGR9 XIRP2_HUMAN          NA
##                                       ProteinName OrganismName GeneName
##                                       <character>        <Rle>    <Rle>
## 1   Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 2   Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 3   Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 4   Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 5   Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## ...                                           ...          ...      ...
## 32  Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 33  Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 34  Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 35  Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
## 36  Xin actin-binding repeat-containing protein 2 Homo sapiens    XIRP2
##              ProteinExistence SequenceVersion Comment  spectrumid
##                         <Rle>           <Rle>   <Rle> <character>
## 1   Evidence at protein level               2      NA   index=124
## 2   Evidence at protein level               2      NA    index=28
## 3   Evidence at protein level               2      NA    index=20
## 4   Evidence at protein level               2      NA    index=21
## 5   Evidence at protein level               2      NA   index=187
## ...                       ...             ...     ...         ...
## 32  Evidence at protein level               2      NA    index=87
## 33  Evidence at protein level               2      NA    index=99
## 34  Evidence at protein level               2      NA     index=9
## 35  Evidence at protein level               2      NA   index=122
## 36  Evidence at protein level               2      NA    index=77
##     acquisitionnum passthreshold      rank calculatedmasstocharge
##          <numeric>     <logical> <integer>            <character>
## 1              124          TRUE         1        715.03076171875
## 2               28          TRUE         1      715.4117431640625
## 3               20          TRUE         1      786.9080810546875
## 4               21          TRUE         1       524.941162109375
## 5              187          TRUE         1      629.3385620117188
## ...            ...           ...       ...                    ...
## 32              87          TRUE         1      720.3527221679688
## 33              99          TRUE         1      618.7781982421875
## 34               9          TRUE         1          1013.51171875
## 35             122          TRUE         1       820.890869140625
## 36              77          TRUE         1      821.9254150390625
##     experimentalmasstocharge chargestate ms.gf.denovoscore ms.gf.evalue
##                  <character>   <integer>         <numeric>    <numeric>
## 1          715.0304565429688           3               146     0.002078
## 2          715.9176635742188           2                12     0.065243
## 3          786.9065551757812           2                35     0.119470
## 4          524.9410400390625           3                49     0.524066
## 5          629.8379516601562           2                53     0.339940
## ...                      ...         ...               ...          ...
## 32         720.3445434570312           2                63        42.38
## 33          619.288818359375           2                52        21.54
## 34        1014.0198364257812           2                17        21.91
## 35         821.4005126953125           2                56        22.61
## 36          821.923095703125           2               102        42.87
##     ms.gf.rawscore ms.gf.specevalue assumeddissociationmethod
##          <numeric>        <numeric>               <character>
## 1               36        1.823e-07                       CID
## 2              -62        5.738e-06                       CID
## 3              -44        1.050e-05                       CID
## 4              -19        4.606e-05                       CID
## 5              -33        2.993e-05                       CID
## ...            ...              ...                       ...
## 32             -90         0.003725                       CID
## 33             -63         0.001895                       CID
## 34            -148         0.001922                       CID
## 35             -97         0.001987                       CID
## 36            -102         0.003766                       CID
##     ctermioncurrentratio explainedioncurrentratio isotopeerror
##              <character>              <character>  <character>
## 1             0.18311569                0.3733551            0
## 2             0.17072089               0.17917182            1
## 3            0.081987105               0.15870962            0
## 4             0.05622067               0.14561509            0
## 5            0.082444645               0.16891155            1
## ...                  ...                      ...          ...
## 32           0.037362944              0.063179314            0
## 33           0.062403493               0.09917885            1
## 34           0.062099736               0.11363181            1
## 35           0.030872315              0.047606274            1
## 36           0.031126814              0.036955375            0
##     meanerrorall meanerrortop7 meanrelerrorall meanrelerrortop7
##      <character>   <character>     <character>      <character>
## 1      57.312866     55.738667       14.757639        46.320156
## 2      340.82858     151.35628        319.4936       126.973404
## 3      48.116673     48.116673       2.6420627        2.6420627
## 4      132.30872     132.30872      0.31386718       0.31386718
## 5      37.566624     37.566624        8.980245         8.980245
## ...          ...           ...             ...              ...
## 32     241.91089     241.91089       148.93636        148.93636
## 33     247.79701     247.79701      -103.02391       -103.02391
## 34     117.01554     117.01554       117.01554        117.01554
## 35     55.865704     55.865704        52.49654         52.49654
## 36     276.21368     276.21368       162.91066        162.91066
##     ms2ioncurrent ntermioncurrentratio nummatchedmainions stdeverrorall
##       <character>          <character>        <character>   <character>
## 1       1984657.5            0.1902394                 10      87.09684
## 2       1164420.1          0.008450932                  8      523.5172
## 3        439678.0            0.0767225                  7      34.59512
## 4       461579.66           0.08939442                  5      93.69792
## 5       1574020.8           0.08646689                  6     36.306866
## ...           ...                  ...                ...           ...
## 32       935006.0           0.02581637                  6     374.24667
## 33      1784755.5           0.03677536                  7     319.73682
## 34      317826.88           0.05153207                  5      76.52627
## 35      1221484.8          0.016733961                  3      68.27473
## 36      829413.75         0.0058285617                  4     213.22662
##     stdeverrortop7 stdevrelerrorall stdevrelerrortop7             pepseq
##        <character>      <character>       <character>        <character>
## 1         92.53631        103.21258          97.59206 QEITQNKSFFSSVKESQR
## 2         161.3231         536.8037         181.13979       LPVPKDVYSKQR
## 3         34.59512        59.203514         59.203514      EQNNDALEKSLRR
## 4         93.69792        162.12587         162.12587      EQNNDALEKSLRR
## 5        36.306866        51.466446         51.466446         SLKESSHRWK
## ...            ...              ...               ...                ...
## 32       374.24667         419.9993          419.9993      TNTSTGLKMAMER
## 33       319.73682        391.17908         391.17908       PESGFAEDSAAR
## 34        76.52627         76.52627          76.52627 QPDAIPGDIEKAIECLEK
## 35        68.27473         70.89802          70.89802     MARYQAAVSRGDCR
## 36       213.22662         308.5769          308.5769     QEITQNKSFFSSVK
##      modified      modification   isdecoy        post         pre
##     <logical>       <character> <logical> <character> <character>
## 1       FALSE                NA     FALSE           D           K
## 2       FALSE                NA     FALSE           N           R
## 3       FALSE                NA     FALSE           L           R
## 4       FALSE                NA     FALSE           L           R
## 5       FALSE                NA     FALSE           E           K
## ...       ...               ...       ...         ...         ...
## 32      FALSE                NA     FALSE           S           K
## 33      FALSE                NA     FALSE           G           K
## 34       TRUE 57.021463735 (15)     FALSE           A           K
## 35       TRUE 57.021463735 (13)     FALSE           S           R
## 36      FALSE                NA     FALSE           E           K
##           end     start    length                  spectrumFile
##     <integer> <integer> <integer>                   <character>
## 1        2760      2743      3374 Thermo_Hela_PRTC_selected.mgf
## 2         318       307      3374 Thermo_Hela_PRTC_selected.mgf
## 3        1870      1858      3374 Thermo_Hela_PRTC_selected.mgf
## 4        1870      1858      3374 Thermo_Hela_PRTC_selected.mgf
## 5        1708      1699      3374 Thermo_Hela_PRTC_selected.mgf
## ...       ...       ...       ...                           ...
## 32       2094      2082      3374 Thermo_Hela_PRTC_selected.mgf
## 33         31        20      3374 Thermo_Hela_PRTC_selected.mgf
## 34       1729      1712      3374 Thermo_Hela_PRTC_selected.mgf
## 35         61        48      3374 Thermo_Hela_PRTC_selected.mgf
## 36       2756      2743      3374 Thermo_Hela_PRTC_selected.mgf
##                    databaseFile
##                     <character>
## 1   HUMAN_2014_06-09prots.fasta
## 2   HUMAN_2014_06-09prots.fasta
## 3   HUMAN_2014_06-09prots.fasta
## 4   HUMAN_2014_06-09prots.fasta
## 5   HUMAN_2014_06-09prots.fasta
## ...                         ...
## 32  HUMAN_2014_06-09prots.fasta
## 33  HUMAN_2014_06-09prots.fasta
## 34  HUMAN_2014_06-09prots.fasta
## 35  HUMAN_2014_06-09prots.fasta
## 36  HUMAN_2014_06-09prots.fasta
## 
## ...
## <5 more elements>
```

Specific proteins can be extracted by index of name using
`[` and proteins and their peptide features can be plotted
with the default plot method.


```r
seqnames(p)
```

```
## [1] "A4UGR9"   "A6H8Y1"   "O43707"   "O75369"   "P00558"   "P02545"  
## [7] "P04075"   "P04075-2" "P60709"
```

```r
plot(p[c(1,9)])
```

<img src="figure/pplot.png" title="plot of chunk pplot" alt="plot of chunk pplot" style="display: block; margin: auto;" />

More details can be found in `?Proteins`. The object generated above
is also directly available as `data(p)`.

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
## [1] grid      parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] Pbase_0.1.5           Gviz_1.9.10           GenomicRanges_1.17.23
##  [4] GenomeInfoDb_1.1.12   IRanges_1.99.22       S4Vectors_0.1.2      
##  [7] Rcpp_0.11.2           BiocGenerics_0.11.3   rmarkdown_0.2.46     
## [10] knitr_1.6            
## 
## loaded via a namespace (and not attached):
##  [1] affy_1.43.3               affyio_1.33.0            
##  [3] AnnotationDbi_1.27.8      BatchJobs_1.3            
##  [5] BBmisc_1.7                Biobase_2.25.0           
##  [7] BiocInstaller_1.15.5      BiocParallel_0.7.7       
##  [9] biomaRt_2.21.1            Biostrings_2.33.12       
## [11] biovizBase_1.13.8         bitops_1.0-6             
## [13] brew_1.0-6                BSgenome_1.33.8          
## [15] checkmate_1.1             cleaver_1.3.7            
## [17] cluster_1.15.2            codetools_0.2-8          
## [19] colorspace_1.2-4          data.table_1.9.2         
## [21] DBI_0.2-7                 dichromat_2.0-0          
## [23] digest_0.6.4              doParallel_1.0.8         
## [25] evaluate_0.5.5            fail_1.2                 
## [27] foreach_1.4.2             formatR_0.10             
## [29] Formula_1.1-1             GenomicAlignments_1.1.20 
## [31] GenomicFeatures_1.17.12   ggplot2_1.0.0            
## [33] gtable_0.1.2              Hmisc_3.14-4             
## [35] htmltools_0.2.4           impute_1.39.0            
## [37] iterators_1.0.7           lattice_0.20-29          
## [39] latticeExtra_0.6-26       limma_3.21.10            
## [41] MALDIquant_1.10           MASS_7.3-33              
## [43] matrixStats_0.10.0        MSnbase_1.13.12          
## [45] munsell_0.4.2             mzID_1.3.2               
## [47] mzR_1.11.7                pcaMethods_1.55.0        
## [49] plyr_1.8.1                preprocessCore_1.27.1    
## [51] proto_0.3-10              Pviz_0.99.0              
## [53] R.methodsS3_1.6.1         RColorBrewer_1.0-5       
## [55] RCurl_1.95-4.1            reshape2_1.4             
## [57] Rsamtools_1.17.31         RSQLite_0.11.4           
## [59] rtracklayer_1.25.13       scales_0.2.4             
## [61] sendmailR_1.1-2           splines_3.2.0            
## [63] stats4_3.2.0              stringr_0.6.2            
## [65] survival_2.37-7           tools_3.2.0              
## [67] VariantAnnotation_1.11.12 vsn_3.33.0               
## [69] XML_3.98-1.1              XVector_0.5.7            
## [71] zlibbioc_1.11.1
```
