context("Functionality to map peptide features to genome using ensembldb")
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
    ## width of GRanges for one feature have to match lenght of its peptide
    ## sequence * 3:
    res_l <- split(res, res$group)
    expect_true(all(unlist(
        lapply(res_l, function(z) {
            sum(width(z)) == nchar(z$pepseq[1]) * 3
        })
    )))

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
    ## width of GRanges for one feature have to match lenght of its peptide
    ## sequence * 3:
    res_l <- split(res, res$group)
    expect_true(all(unlist(lapply(res_l, function(z) {
        sum(width(z)) == nchar(z$pepseq[1]) * 3
    }))))

    ## - strand:
    ## transcript with 3 exons:
    ## ex1: 60-50, ex2: 40-30, ex3: 10-20
    exs <- GRanges(seqnames = "test",
                   ranges = IRanges(start = c(50, 30, 10),
                                    end = c(60, 40, 20)), strand = "-")
    mcols(exs) <- DataFrame(gene = "a", symbol = "b", transcript = "c")
    ## protein with 3 peptide features: first spanning exons 1 and 2, second
    ## within exon 2 and third spanning exon 2 and 3.
    prot_seq <- paste0(rep("Y", sum(width(exs))/3), collapse = "")
    aa <- AAStringSet(x = prot_seq)
    names(aa) <- "test_tx"
    mcols(aa) <- DataFrame(gene = "test_gene")
    ## Features:
    ## pf_1: 3-5, pf_2: 6-7, pf_3: 7-10
    prng <- IRangesList(IRanges(start = c(3, 6, 7), end = c(5, 7, 10)))
    names(prng[[1]]) <- c("pf_1", "pf_2", "pf_3")
    mcols(prng[[1]]) <- DataFrame(source = rep("manual", 3))
    prt <- new("Proteins", aa = aa, pranges = prng)
    ## The tricky thing with the negative strand is that we have to swith
    ## start and end coordinates of the exons to calculate from positions within
    ## the cDNA to DNA.
    res <- Pbase:::.mapToGenome2(prt, exs)
    ## res <- Pbase:::.mapToGenome(prt, sort(exs))
    ## We expect:
    ## pf_1: 37-40, 50-54
    ## pf_2: 31-36
    ## pf_3: 13-20, 30-33
    ## The order will be pf_3, pf_2, pf_1
    expect_equal(start(res), c(13, 30, 31, 37, 50))
    expect_equal(end(res), c(20, 33, 36, 40, 54))
    ## Again check that rev-splicing was correct
    res_l <- split(res, res$group)
    expect_true(all(unlist(lapply(res_l, function(z) {
        sum(width(z)) == nchar(z$pepseq[1]) * 3
    }))))
})
