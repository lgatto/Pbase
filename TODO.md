
- UniProt fasta header parsing

- A `seqnames` method to access the protein sequence names? Maybe
  `seqids`? Let's see what the UniProt header parser does. Note that
  the seqnames generic comes from GenomicRanges... we could possibly
  ask to move it to `BiocGenerics`? Alternatively, we could use
  `names`, which is used in `Biostrings`. Difference between
  `seqnames` and `names` should be clarified.


