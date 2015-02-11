---
title: "Mapping between Ensembl and UCSC"
output:
  BiocStyle::html_document:
    toc: true
---


<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Ensembl and UCSC mapping}
%\VignettePackage{Pbase}
-->

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").style.marginTop = "0";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>
<style type="text/css" scoped>
body, td {
   font-family: sans-serif;
   background-color: white;
   font-size: 13px;
}

body {
  max-width: 800px;
  margin: 0 auto;
  padding: 1em 1em 2em;
  line-height: 20px;
}

/* element spacing */

p, pre { 
  margin: 0em 0em 1em;
}

/* center images and tables */
img, table {
  margin: 0em auto 1em;
}

p {
  text-align: justify;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, monospace;
}

h1, h2, h3, h4, h5, h6 { 
  font-family: Helvetica, Arial, sans-serif;
  margin: 1.2em 0em 0.6em 0em;
  font-weight: bold;
}

h1 {
  font-size: 250%;
  font-weight: normal;
  color: #87b13f;
  line-height: 1.1em;
}

h2 {
  font-size: 160%;
  font-weight: normal;
  line-height: 1.4em;
  border-bottom: 1px #1a81c2 solid;
}

h3 {
  font-size: 130%;  
}

h2, h3 {
  color: #1a81c2;
}

h4, h5, h6 {
  font-size:115%;
} /* not expecting to dive deeper than four levels on a single page */

/* links are simply blue, hovering slightly less blue */
a { color: #1a81c2; }
a:active { outline: none; }
a:visited { color: #1a81c2; }
a:hover { color: #4c94c2; }

pre, img {
  max-width: 100%;
  display: block;
}

pre {
  border: 0px none;
  background-color: #F8F8F8;
  white-space: pre;
  overflow-x: auto;
}

pre code {
  border: 1px #aaa dashed;
  background-color: white;
  display: block;
  padding: 1em;  
  color: #111;
  overflow-x: inherit;
}

/* markdown v1 */
pre code[class] {
  background-color: inherit;
}

/* markdown v2 */
pre[class] code {
  background-color: inherit;
}

/* formatting of inline code */
code { 
  background-color: transparent;
  color: #87b13f;
  font-size: 92%;
}

/* formatting of tables */

table, td, th {
  border: none;
  padding: 0 0.5em;
}

/* alternating row colors */
tbody tr:nth-child(odd) td {
  background-color: #F8F8F8;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

@media print {
   * {
      background: transparent !important;
      color: black !important;
      filter:none !important;
      -ms-filter: none !important;
   }

   body {
      font-size:12pt;
      max-width:100%;
   }

   a, a:visited {
      text-decoration: underline;
   }

   hr {
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote {
      padding-right: 1em;
      page-break-inside: avoid;
   }

   tr, img {
      page-break-inside: avoid;
   }

   img {
      max-width: 100% !important;
   }

   @page :left {
      margin: 15mm 20mm 15mm 10mm;
   }

   @page :right {
      margin: 15mm 10mm 15mm 20mm;
   }

   p, h2, h3 {
      orphans: 3; widows: 3;
   }

   h2, h3 {
      page-break-after: avoid;
   }
}
</style>

<hr />

Author: [Laurent Gatto](http://cpu.sysbiol.cam.ac.uk/)

Compilation date: 2015-02-11

<hr />


This vignette described how to convert coordinates between Ensembl
release 78 (based on GRCh38) and the UCSC (based on GRCh37). We will
use transcript `ENST00000373316` as working example.


```r
tr <- "ENST00000373316"
```

## Ensembl (GRCh38)

Here is use *[Gviz](http://bioconductor.org/packages/release/bioc/html/Gviz.html)* to query the latest Ensembl biomart
and extract the transcript of interest.


```r
suppressMessages(library("Gviz"))
suppressMessages(library("biomaRt"))

ens <- useMart("ensembl", "hsapiens_gene_ensembl")
etr <- BiomartGeneRegionTrack(biomart = ens,
                              transcript = tr,
                              genome = "hg19")
etr <- split(etr, transcript(etr))
etr <- etr[[tr]]
etr <- ranges(etr)
```

Note the starting position of the transcript is 78104174.

## UCSC (GRCh37)

Below, I repeat the same operation without using my own ens Mart
instance. As far as I understand, Gviz queries the UCSC genome
reference by default. 


```r
utr <- BiomartGeneRegionTrack(transcript = tr,
                              genome = "hg19")
utr <- split(utr, transcript(utr))
utr <- utr[[tr]]
utr <- ranges(utr)
```

Note the starting position of the transcript is 77359671.

These differences seem to stem from different genome builds. **Ensembl
release 78** uses **GRCh38**, while **UCSC** uses **GRCh37**. Indeed,
*[Gviz](http://bioconductor.org/packages/release/bioc/html/Gviz.html)* sets the Ensembl biomart server to `Feb.2014`
`GRCh37.p13`.

## Coordinates conversion

We will use the coordinate mapping infrastructure described in the
[January 2015 Bioconductor Newletter](http://www.bioconductor.org/help/newsletters/2015_January/#coordinate-mapping)
and the
[Changing genomic coordinate systems with rtracklayer::liftOver](http://bioconductor.org/help/workflows/liftOver/)
workflow.

First, we query *[AnnotationHub](http://bioconductor.org/packages/release/bioc/html/AnnotationHub.html)* for a chain file to
perform the operation we want.


```r
library("AnnotationHub")
hub <- AnnotationHub()
query(hub, 'hg19ToHg38')
```

```
## class: AnnotationHub 
## hub: https://annotationhub.bioconductor.org 
## cache: /home/lg390/.AnnotationHub 
## -------------------------------------------------------------- 
## Some common metadata fields represented here (e.g.): 
## Title:  hg19ToHg38.over.chain.gz 
## dataprovider:  hgdownload.cse.ucsc.edu 
## species:  Homo sapiens 
## taxonomyid:  9606 
## genome:  hg19 
## description:  UCSC liftOver chain file from hg19 to hg38 
## tags:  liftOver, chain, UCSC, genome, homology 
## rdataclass:  ChainFile 
## recipe:  NA 
## sourceurl:  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz 
## To see the full range for any of these,  please use the '$' operator and hit tab from this object.
```

```r
chain <- query(hub, 'hg19ToHg38')[[1]]
```

The `liftOver` function from the *[rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html)* package
will use the chain and translate the coordinates of a `GRanges` object
into a new `GRangesList` object.


```r
library("rtracklayer")

res <- liftOver(utr, chain)
res <- unlist(res)

## set annotation
genome(res) <- "hg19"
names(res) <- NULL

identical(res, etr)
```

```
## [1] TRUE
```

## Session information


```r
sessionInfo()
```

```
## R Under development (unstable) (2015-01-22 r67580)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.1 LTS
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
##  [1] grid      parallel  stats4    methods   stats     graphics  grDevices
##  [8] utils     datasets  base     
## 
## other attached packages:
##  [1] rtracklayer_1.27.7    AnnotationHub_1.99.45 XVector_0.7.4        
##  [4] biomaRt_2.23.5        Gviz_1.11.11          GenomicRanges_1.19.36
##  [7] GenomeInfoDb_1.3.12   IRanges_2.1.38        S4Vectors_0.5.19     
## [10] BiocGenerics_0.13.4   BiocStyle_1.5.3       rmarkdown_0.5.1      
## [13] knitr_1.9            
## 
## loaded via a namespace (and not attached):
##  [1] acepack_1.3-3.3              AnnotationDbi_1.29.17       
##  [3] base64enc_0.1-2              BatchJobs_1.5               
##  [5] BBmisc_1.9                   Biobase_2.27.1              
##  [7] BiocInstaller_1.17.5         BiocParallel_1.1.13         
##  [9] Biostrings_2.35.7            biovizBase_1.15.2           
## [11] bitops_1.0-6                 brew_1.0-6                  
## [13] BSgenome_1.35.16             checkmate_1.5.1             
## [15] cluster_2.0.1                codetools_0.2-10            
## [17] colorspace_1.2-4             DBI_0.3.1                   
## [19] dichromat_2.0-0              digest_0.6.8                
## [21] evaluate_0.5.5               fail_1.2                    
## [23] foreach_1.4.2                foreign_0.8-62              
## [25] formatR_1.0                  Formula_1.2-0               
## [27] GenomicAlignments_1.3.27     GenomicFeatures_1.19.18     
## [29] Hmisc_3.14-6                 htmltools_0.2.6             
## [31] httpuv_1.3.2                 httr_0.6.1                  
## [33] interactiveDisplayBase_1.5.1 iterators_1.0.7             
## [35] lattice_0.20-29              latticeExtra_0.6-26         
## [37] matrixStats_0.13.1           mime_0.2                    
## [39] munsell_0.4.2                nnet_7.3-8                  
## [41] plyr_1.8.1                   R6_2.0.1                    
## [43] RColorBrewer_1.1-2           Rcpp_0.11.4                 
## [45] RCurl_1.95-4.5               RJSONIO_1.3-0               
## [47] R.methodsS3_1.6.1            rpart_4.1-8                 
## [49] Rsamtools_1.19.27            RSQLite_1.0.0               
## [51] scales_0.2.4                 sendmailR_1.2-1             
## [53] shiny_0.11                   splines_3.2.0               
## [55] stringr_0.6.2                survival_2.37-7             
## [57] tools_3.2.0                  VariantAnnotation_1.13.28   
## [59] XML_3.98-1.1                 xtable_1.7-4                
## [61] zlibbioc_1.13.0
```
