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

    ## issue #29
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


test_that("mapToGenome,Proteins,EnsDb", {
    library(EnsDb.Hsapiens.v86)
    edb <- EnsDb.Hsapiens.v86
    ## Let's start with our famous ZBTB16 example:
    zbtb16 <- Proteins(edb, filter = GenenameFilter("ZBTB16"))
    zbtb16_cds <- cdsBy(edb,
                        filter = TxidFilter(acols(zbtb16)$tx_id),
                        columns = c("tx_id", "protein_id"))
    zbtb16_cds <- zbtb16_cds[acols(zbtb16)$tx_id]

    ## Do the mapping using EnsDb, using the protein_id
    res <- mapToGenome(zbtb16, edb)
    ## Check if we get what we expect:
    for (i in 1:length(res)) {
        ## tx_id, protein_id from Proteins and result have to match
        expect_equal(unique(res[[i]]$tx_id), acols(zbtb16[i])$tx_id)
        expect_equal(unique(res[[i]]$accession), seqnames(zbtb16)[i])
        expect_equal(res[[i]],
                    mapToGenome(zbtb16[i], zbtb16_cds[i])[[1]])
    }
    ## Use tx_id as identifiers.
    res_2 <- mapToGenome(zbtb16, edb, id = "tx_id", idType = "tx_id")
    expect_equal(res, res_2)

    ## Use uniprot ID.
    uniprts <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                        columns = "uniprot_id", return.type = "data.frame")
    uniprts <- unique(uniprts$uniprot_id)
    prts <- Proteins(edb, filter = UniprotidFilter(uniprts))
    ## FIX this does not return the uniprot ids.

    ## Check errors:
    expect_error(mapToGenome(zbtb16, edb, id = "other"))
    expect_error(mapToGenome(zbtb16, edb, idType = "exon"))
})


benchmark_pmapToGenome <- function() {
    ## issue #20
    ## Simple benchmark comparing the efficiency of the pmapToGenome method
    ## using sapply and a for loop with the pmapToGenom2 that uses bpmapply.
    library(microbenchmark)
    library(BiocParallel)

    ## Fetch a bunch of proteins.
    gnf <- GenenameFilter(c("ZBTB16", "BCL2", "BCL2L11"))
    cdss <- cdsBy(edb, by = "tx", columns = c("tx_id", "protein_id"),
                  filter = gnf)
    prts <- Proteins(edb, filter = gnf)
    cdss <- cdss[acols(prts)$tx_id]

    ## (1) Without parallel processing
    ## FIX THE ERRORS:
    register(SerialParam())
    res <- pmapToGenome(prts, cdss)
    res_2 <- Pbase:::pmapToGenome2(prts, cdss)
    expect_equal(res, res_2)
    microbenchmark(pmapToGenome(prts, cdss),
                   Pbase:::pmapToGenome2(prts, cdss),
                   times = 20)
    ## Unit: seconds
    ##                               expr      min       lq     mean   median       uq
    ##           pmapToGenome(prts, cdss) 1.334235 1.371089 1.537094 1.442745 1.610697
    ##  Pbase:::pmapToGenome2(prts, cdss) 1.624151 1.644120 1.797967 1.672039 1.744216
    ##       max neval cld
    ##  2.272202    20  a
    ##  3.827477    20   b

    ## With 2 CPUs.
    register(MulticoreParam(workers = 2))
    microbenchmark(pmapToGenome(prts, cdss),
                   Pbase:::pmapToGenome2(prts, cdss),
                   times = 20)
    ## Unit: seconds
    ##                               expr      min       lq     mean   median       uq
    ##           pmapToGenome(prts, cdss) 1.348400 1.370049 1.520549 1.453795 1.502679
    ##  Pbase:::pmapToGenome2(prts, cdss) 4.108021 4.148563 4.252193 4.163558 4.213703
    ##       max neval cld
    ##  2.218503    20  a
    ##  5.573078    20   b
}

benchmark_dot_mapToGenome <- function() {
    ## issue #30
    ## Benchmark comparing the speed of the .mapToGenome and .mapToGenome2
    ## functions.
    ## Need an object that is compatible with both, .mapToGenoma and
    ## .mapToGenom2, i.e. has mcols 'gene', 'transcript' and 'symbol' in the
    ## GRanges
    zbtb16 <- Proteins(edb, filter = GenenameFilter("ZBTB16"))
    ## NOTE: the first query is much faster!
    zbtb16_cds <- cdsBy(edb,
                        filter = TxidFilter(acols(zbtb16)$tx_id),
                        columns = c("tx_id", "protein_id"))
    prot <- zbtb16[1]
    gen <- zbtb16_cds[acols(prot)$tx_id]
    gen <- unlist(gen)
    gen$gene <- "ZBTB16"
    gen$symbol <- "ZBTB16"
    gen$transcript <- gen$tx_id

    res <- Pbase:::.mapToGenome(prot, gen)
    res_2 <- Pbase:::.mapToGenome2(prot, gen)
    ## .mapToGenome2 returns more mcols and names, fix that.
    res <- unlist(res)
    res_2 <- unname(unlist(res_2))
    mcols(res_2) <- mcols(res_2)[, colnames(mcols(res))]
    ## Drop column group as this is different.
    mcols(res) <- mcols(res)[, colnames(mcols(res)) != "group"]
    mcols(res_2) <- mcols(res_2)[, colnames(mcols(res_2)) != "group"]
    expect_equal(res, res_2)
    ## OK

    ## Benchmark:
    library(microbenchmark)
    microbenchmark(Pbase:::.mapToGenome(prot, gen),
                   Pbase:::.mapToGenome2(prot, gen), times = 20)
    ## Unit: milliseconds
    ##                              expr       min        lq      mean    median
    ##   Pbase:::.mapToGenome(prot, gen) 918.21352 927.42580 972.29711 940.74852
    ##  Pbase:::.mapToGenome2(prot, gen)  50.62116  53.04321  55.59854  55.08398
    ##         uq        max neval cld
    ##  952.77226 1583.10086    20   b
    ##   56.63087   67.59263    20  a
}


## TODO:
## + Add more unit tests for the mapToGenome,Proteins,EnsDb (uniprot, tx_id etc)
## + issue: Proteins with uniprot_id: problem is the 1:n mapping from protein_id
##   to uniprot_id in Ensembl. Solution is to return a redundant list of proteins
##   i.e. replicate the protein_id, one for each uniprot_id.
