contect("Functionality to map peptide features to genome using ensembldb")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

test_that(".mapToGenome2 internal function", {
    library(ensembldb)
    library(Biostrings)
    ## Real life example.
    zbtb16 <- Proteins(edb, filter = GenenameFilter("ZBTB16"))
    zbtb16_cds <- cdsBy(edb,
                        filter = TxidFilter(acols(zbtb16)$tx_id))
    zbtb16_cds <- zbtb16_cds[acols(zbtb16)$tx_id]
    res <- Pbase:::.mapToGenome2(zbtb16[1], grObj = zbtb16_cds[[1]])
    ## Check the mcols.
    expect_true(all(res$tx_id == acols(zbtb16)$tx_id[1]))
    ## Check that names are correct and match the sequence:
    pfeat <- pfeatures(zbtb16)[[1]]
    ## All sequences have to be present in res:
    expect_true(all(as.character(pfeat) %in% res$pepseq))
    ## The names have to match the peptide sequences.
    pfeat_mat <- cbind(name = names(pfeat), seq = as.character(pfeat))
    ## There can be a 1:n mapping between sequence and ID
    pfeat_list <- split(pfeat_mat[, 1], pfeat_mat[, 2])
    for (i in 1:length(res)) {
        expect_true(names(res)[i] %in% pfeat_list[[res$pepseq[i]]])
    }

    ## Changing mcols of the GRanges and check that we get what we expect.
    tmp <- unlist(zbtb16_cds)
    ## Empty mcols
    mcols(tmp) <- NULL
    tmp <- split(tmp, unlist(zbtb16_cds)$tx_id)
    res <- Pbase:::.mapToGenome2(zbtb16[1], grObj = tmp[[1]])
    ## Expect only columns pepseq, accession, exonJunctions, group in mcols
    expect_equal(colnames(mcols(res)),
                 c("pepseq", "accession", "exonJunctions", "group"))
    ## Now add one columns "symbol"
    tmp <- unlist(tmp)
    mcols(tmp) <- DataFrame(symbol = rep("ZBTB16", length(tmp)))
    tmp <- split(tmp, unlist(zbtb16_cds)$tx_id)
    res <- Pbase:::.mapToGenome2(zbtb16[1], grObj = tmp[[1]])
    ## Expect column symbol also in res:
    expect_true(all(res$symbol == "ZBTB16"))

    ## Artificial examples to check that mapping works for forward and reverse
    ## strand.
    ## + strand:
    ## transcript with 3 exons:
    ## 11-20, 25-40, 50-80
    exs <- GRanges(seqnames = "j",
                   ranges = IRanges(start = c(11, 25, 50),
                                    end = c(20, 40, 80)), strand = "+")
    ## protein with 2 peptide features, first spanning exons 1 and 2, second in
    ## exon 3.
    prot_seq <- paste0(rep("Y", sum(width(exs))/3), collapse = "")
    aa <- AAStringSet(x = prot_seq)
    names(aa) <- "test_tx"
    mcols(aa) <- DataFrame(tx_id = "test_tx")
    ## features: 2-8 and 15-18
    prng <- IRangesList(IRanges(start = c(2, 15), end = c(8, 18)))
    names(prng[[1]]) <- c("pf_1", "pf_2")
    mcols(prng[[1]]) <- DataFrame(source = c("manual", "manual"))
    prt <- new("Proteins", aa = aa, pranges = prng)
    ## What do we expect:
    ## first peptide feature: 14-20, 25-38
    ## second peptide feature: 66-77
    res <- Pbase:::.mapToGenome2(prt, exs)
    expect_equal(start(res), c(14, 25, 66))
    expect_equal(end(res), c(20, 38, 77))

    ## - strand:
})
