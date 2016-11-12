context("Proteins methods and functions related to ensembldb")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

test_that("Proteins,EnsDb,missing constructor", {
    library(ensembldb)
    library(RSQLite)
    ## Load Proteins with protein domains.
    prots <- Proteins(edb, filter = GenenameFilter("ZBTB16"))
    ## Check that we've got all proteins.
    txs <- transcripts(edb, filter = GenenameFilter("ZBTB16"))
    res <- dbGetQuery(dbconn(edb),
                      paste0("select * from protein where tx_id in (",
                             paste0("'", txs$tx_id, "'", collapse = ", "),")"))
    ## protein_id
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    ## tx_id
    expect_identical(sort(acols(prots)$tx_id), sort(res$tx_id))
    ## sequence
    expect_identical(unname(as.character(aa(prots)[res$protein_id])),
                     res$protein_sequence)
    ## Check that we have all protein domains and that they are correctly
    ## assigned to their respective protein.
    expect_identical(seqnames(prots), names(pranges(prots)))
    ## All nrow pcols except 1 should be > 0
    expect_true(sum(unlist(lapply(pcols(prots), nrow)) == 0) == 1)
    ## Check that the mcols also matches:
    expect_true(length(pcols(prots)) == length(unique(res$protein_id)))
    ## Check that we've got all protein domains
    pd <- dbGetQuery(dbconn(edb),
                     paste0("select * from protein_domain where protein_id in (",
                            paste0("'", res$protein_id,"'", collapse = ", "),
                            ")"))
    pdL <- split(pd, pd$protein_id)
    tmp <- lapply(pdL, function(z) {
        tmp <- pranges(prots)[[z$protein_id[1]]]
        expect_identical(start(tmp), z$prot_dom_start)
        expect_identical(mcols(tmp)$interpro_accession, z$interpro_accession)
    })
    ## Without protein domains.
    prots <- Proteins(file = edb, filter = GenenameFilter("ZBTB16"),
                      loadProteinDomains = FALSE)
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    expect_identical(seqnames(prots), names(pranges(prots)))
    ## Lenght (nrow) of pranges/pcols is 0
    expect_true(all(unlist(lapply(pcols(prots), nrow)) == 0))
    expect_true(all(unlist(lengths(pranges(prots))) == 0))

    ## Other gene: BCL2L11
    prots <- Proteins(edb, filter = GenenameFilter("BCL2L11"), fetchLRG = TRUE)
    ## Check that we've got all proteins.
    txs <- transcripts(edb, filter = GenenameFilter("BCL2L11"))
    res <- dbGetQuery(dbconn(edb),
                      paste0("select * from protein where tx_id in (",
                             paste0("'", txs$tx_id, "'", collapse = ", "),")"))
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    expect_identical(sort(acols(prots)$tx_id), sort(res$tx_id))
    expect_identical(unname(as.character(aa(prots)[res$protein_id])),
                     res$protein_sequence)
    expect_identical(seqnames(prots), names(pranges(prots)))
    expect_identical(lengths(pranges(prots)),
                     unlist(lapply(pcols(prots), nrow)))
    pd <- dbGetQuery(dbconn(edb),
                     paste0("select * from protein_domain where protein_id in (",
                            paste0("'", res$protein_id,"'", collapse = ", "),
                            ")"))
    pdL <- split(pd, pd$protein_id)
    tmp <- lapply(pdL, function(z) {
        tmp <- pranges(prots)[[z$protein_id[1]]]
        expect_identical(start(tmp), z$prot_dom_start)
        expect_identical(mcols(tmp)$interpro_accession, z$interpro_accession)
    })
    ## Without protein domains.
    prots <- Proteins(file = edb, filter = GenenameFilter("BCL2L11"),
                      loadProteinDomains = FALSE, fetchLRG = TRUE)
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    expect_identical(seqnames(prots), names(pranges(prots)))
    ## Lenght (nrow) of pranges/pcols is 0
    expect_true(all(unlist(lapply(pcols(prots), nrow)) == 0))
    expect_true(all(unlist(lengths(pranges(prots))) == 0))
})

test_that("Proteins,EnsDb,missing protein_id n:1 uniprot_id mapping", {
    dat <- proteins(edb, filter = GenenameFilter("ZBTB16"),
                    columns = c("tx_id", "protein_id", "uniprot_id"),
                    return.type = "data.frame")
    uniprts <- unique(dat$uniprot_id)
    prts <- Proteins(edb, filter = UniprotidFilter(uniprts))

    ## Check content; ordering should be the same.
    expect_equal(seqnames(prts), dat$protein_id)
    expect_equal(acols(prts)$tx_id, dat$tx_id)
    expect_equal(acols(prts)$uniprot_id, dat$uniprot_id)
    ## Protein domains.
    for (i in 1:length(acols)) {
        res <- proteins(edb, filter = ProteinidFilter(seqnames(prts)[i]),
                        columns = listColumns(edb, "protein_domain"))
        expect_equal(start(pranges(prts)[[i]]), res$prot_dom_start)
        expect_equal(end(pranges(prts)[[i]]), res$prot_dom_end)
        expect_equal(mcols(pranges(prts)[[i]])$interpro_accession,
                     res$interpro_accession)
    }
})

dontrun_test_multi_unitprot <- function() {
    library(EnsDb.Hsapiens.v86)
    edb <- EnsDb.Hsapiens.v86
    gn <- "XIRP2" ## that's the first gene from the data(p)
    prts <- proteins(edb, filter = GenenameFilter(gn),
                     columns = c("tx_id","uniprot_id"))
    ## OK, here I've got multiple proteins for one Uniprot.
    res <- Proteins(edb, filter = GenenameFilter(gn),
                    columns = c("tx_id", "gene_name", "uniprot_id"))
    res_2 <- Proteins(edb, filter = UniprotidFilter("A4UGR9"),
                      columns = c("tx_id", "gene_name"))
}

dontrun_test_all <- function() {
    ## Get the mapping between all protein_id and tx_id
    library(RSQLite)
    maps <- dbGetQuery(dbconn(edb), "select protein_id,tx_id from protein;")
    ## Have 103,725 entries, protein_id: 103,722, tx_id: 103,725.
    head(sort(table(maps$protein_id), decreasing = TRUE))
    ## OK, the LRG have multi-mappings! LRG_321p8:3, LRG_321p1:2
    ## Fetching all proteins EXCEPT LRG proteins from the database:
    system.time(
        all_prots <- Proteins(edb, filter = TxidFilter("ENS%",
                                                       condition = "like"),
                              loadProteinDomains = FALSE)
    ) ## 20 sec
    system.time(
        all_prots <- Proteins(edb, fetchLRG = TRUE,
                              loadProteinDomains = FALSE)
    )
    system.time(
        all_prots <- Proteins(edb, filter = TxidFilter("ENS%",
                                                       condition = "like"),
                              loadProteinDomains = TRUE)
    ) ## 93 sec
    ## Seems like we have to use the TxidFilter by default to avoid fetching LRG
    ## genes.
}

dontrun_IRangesList_unique <- function() {
    ## Simple test just to evaluate whether IRanges and IRangesList can have
    ## duplicated names.

    ## GRanges with non-unique names.
    ir <- IRanges(start = c(4, 5, 10), end = c(12, 35, 34))
    library(GenomicFeatures)
    gr <- GRanges(seqnames = rep(1, 3), ranges = ir)
    mcols(gr) <- DataFrame(tx_id = c("tx_1", "tx_2", "tx_1"))
    ## Let's see if names can be non-unique:
    names(gr) <- c("a", "b", "a")
    ## that's OK.
    grL <- split(gr, gr$tx_id)
    ## Let's see if names of GRangesList can be non-unique
    names(grL) <- c("tx_3", "tx_3")
    ## that's OK too
    names(grL)
}
