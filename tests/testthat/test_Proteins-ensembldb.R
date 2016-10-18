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
    prots <- Proteins(edb, filter = GenenameFilter("BCL2L11"))
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
                      loadProteinDomains = FALSE)
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    expect_identical(seqnames(prots), names(pranges(prots)))
    ## Lenght (nrow) of pranges/pcols is 0
    expect_true(all(unlist(lapply(pcols(prots), nrow)) == 0))
    expect_true(all(unlist(lengths(pranges(prots))) == 0))
})
