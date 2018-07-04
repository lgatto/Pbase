context("Proteins methods and functions related to ensembldb")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

test_that("Proteins,EnsDb,missing constructor", {
    library(ensembldb)
    library(RSQLite)
    ## Load Proteins with protein domains.
    prots <- Proteins(edb, filter = ~ gene_name == "ZBTB16")
    ## Check that we've got all proteins.
    txs <- transcripts(edb, filter = GeneNameFilter("ZBTB16"))
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
    tmp <- proteins(edb, filter = ProteinIdFilter(seqnames(prots)),
                    columns = "tx_id")
    rownames(tmp) <- tmp$protein_id
    expect_identical(tmp[seqnames(prots), "tx_id"], mcols(prots@aa)$tx_id)
    expect_identical(seqnames(prots), names(pcols(prots)[, 1]))
    ## All nrow pcols except 1 should be > 0
    expect_true(sum(lengths(pcols(prots)[, 1]) == 0) == 1)
    ## Check that the mcols also matches:
    expect_true(nrow(pcols(prots)) == length(unique(res$protein_id)))
    ## Check that we've got all protein domains
    pd <- dbGetQuery(dbconn(edb),
                     paste0("select * from protein_domain where protein_id in (",
                            paste0("'", res$protein_id,"'", collapse = ", "),
                            ")"))
    pdL <- split(pd, pd$protein_id)[seqnames(prots)]
    for (i in 1:length(pdL)) {
        tmp <- pranges(prots)[i, ][[1]]
        if (length(tmp) == 0) {
            expect_true(length(pdL[[i]]) == 0)
        } else {
            expect_identical(start(tmp), pdL[[i]]$prot_dom_start)
            expect_identical(mcols(tmp)$interpro_accession,
                             pdL[[i]]$interpro_accession)
        }
    }
    ## Without protein domains.
    prots <- Proteins(file = edb, filter = GeneNameFilter("ZBTB16"),
                      loadProteinDomains = FALSE)
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    expect_identical(length(prots), nrow(pranges(prots)))
    ## Lenght (nrow) of pranges/pcols is 0
    expect_true(all(unlist(lapply(pcols(prots), nrow)) == 0))
    expect_true(all(unlist(lengths(pranges(prots))) == 0))

    ## Other gene: BCL2L11
    prots <- Proteins(edb, filter = GeneNameFilter("BCL2L11"), fetchLRG = TRUE)
    ## Check that we've got all proteins.
    txs <- transcripts(edb, filter = GeneNameFilter("BCL2L11"))
    res <- dbGetQuery(dbconn(edb),
                      paste0("select * from protein where tx_id in (",
                             paste0("'", txs$tx_id, "'", collapse = ", "),")"))
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    expect_identical(sort(acols(prots)$tx_id), sort(res$tx_id))
    expect_identical(unname(as.character(aa(prots)[res$protein_id])),
                     res$protein_sequence)
    rownames(res) <- res$protein_id
    expect_identical(res[seqnames(prots), "tx_id"], mcols(prots@aa)$tx_id)
    expect_identical(seqnames(prots), names(pcols(prots)[, 1]))
    
    pd <- dbGetQuery(dbconn(edb),
                     paste0("select * from protein_domain where protein_id in (",
                            paste0("'", res$protein_id,"'", collapse = ", "),
                            ")"))
    pdL <- split(pd, pd$protein_id)
    tmp <- lapply(pdL, function(z) {
        tmp <- pranges(prots)[, 1][[z$protein_id[1]]]
        expect_identical(start(tmp), z$prot_dom_start)
        expect_identical(mcols(tmp)$interpro_accession, z$interpro_accession)
    })
    ## Without protein domains.
    prots <- Proteins(file = edb, filter = GeneNameFilter("BCL2L11"),
                      loadProteinDomains = FALSE, fetchLRG = TRUE)
    expect_identical(sort(seqnames(prots)), sort(res$protein_id))
    expect_identical(length(prots), nrow(pranges(prots)))
})

test_that("Proteins,EnsDb,missing protein_id n:1 uniprot_id mapping", {
    dat <- proteins(edb, filter = GeneNameFilter("ZBTB16"),
                    columns = c("tx_id", "protein_id", "uniprot_id"),
                    return.type = "data.frame")
    uniprts <- unique(dat$uniprot_id)
    prts <- Proteins(edb, filter = UniprotFilter(uniprts))

    ## Check content; ordering should be the same.
    expect_equal(seqnames(prts), dat$protein_id)
    expect_equal(acols(prts)$tx_id, dat$tx_id)
    expect_equal(acols(prts)$uniprot_id, dat$uniprot_id)
    ## Protein domains.
    for (i in 1:length(acols)) {
        res <- proteins(edb, filter = ProteinIdFilter(seqnames(prts)[i]),
                        columns = listColumns(edb, "protein_domain"))
        expect_equal(start(pranges(prts)[, 1][[i]]), res$prot_dom_start)
        expect_equal(end(pranges(prts)[, 1][[i]]), res$prot_dom_end)
        expect_equal(mcols(pranges(prts)[, 1][[i]])$interpro_accession,
                     res$interpro_accession)
    }
    idmap <- unique(cbind(names(prts), acols(prts)$uniprot_id))
    ## We have an n:m mapping:
    expect_true(length(unique(idmap[, 1])) < nrow(idmap))
    expect_true(length(unique(idmap[, 2])) < nrow(idmap))

    ## We can use a UniprotmappingFilter to restrict to good quality maps.
    prts_2 <- Proteins(edb, filter = AnnotationFilterList(
                                GeneNameFilter("ZBTB16"),
                                UniprotMappingTypeFilter("DIRECT")),
                       columns = c("uniprot_id", "protein_id"))
    ## That way we reduce it to a n:1 mapping between Ensembl protein ID and
    ## Uniprot ID.
    idmap <- unique(cbind(names(prts_2), acols(prts_2)$uniprot_id))
    expect_true(length(unique(idmap[, 1])) == nrow(idmap))
    expect_true(length(unique(idmap[, 2])) < nrow(idmap))
})


dontrun_test_multi_unitprot <- function() {
    library(EnsDb.Hsapiens.v86)
    edb <- EnsDb.Hsapiens.v86
    gn <- "XIRP2" ## that's the first gene from the data(p)
    prts <- proteins(edb, filter = GeneNameFilter(gn),
                     columns = c("tx_id","uniprot_id"))
    ## OK, here I've got multiple proteins for one Uniprot.
    res <- Proteins(edb, filter = GeneNameFilter(gn),
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
        all_prots <- Proteins(edb, filter = TxIdFilter("ENS%",
                                                       condition = "like"),
                              loadProteinDomains = FALSE)
    ) ## 20 sec
    system.time(
        all_prots <- Proteins(edb, fetchLRG = TRUE,
                              loadProteinDomains = FALSE)
    )
    system.time(
        all_prots <- Proteins(edb, filter = TxIdFilter("ENS%",
                                                       condition = "like"),
                              loadProteinDomains = TRUE)
    ) ## 93 sec
    ## Seems like we have to use the TxIdFilter by default to avoid fetching LRG
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
