setMethod("Proteins", c("missing", "missing"),
          function(file, uniprotIds, ...){
              aa <- new("AAStringSet")
              aa@elementMetadata <- DataFrame()
              new("Proteins",
                  metadata = list(created = date()),
                  aa = aa,
                  pranges = IRangesList())
          })

setMethod("Proteins",
          signature(file = "character", uniprotIds = "missing"),
          function(file, uniprotIds, ...) {
              p <- .ProteinsFromFasta(filenames = file, ...)
              names(p@aa) <- acols(p)$AccessionNumber
              p
          })

setMethod("Proteins",
          signature(file = "missing", uniprotIds = "character"),
          function(file, uniprotIds, ...) {
              .toBeImplemented()
          })

############################################################
## Proteins for EnsDb, character, use uniprotIds in an UniprotidFilter.
## Proteins for EnsDb, missing, ... -> check if we've got filter specified,
##  otherwise just return everything.
## Add also an optional argument loadDomains = TRUE -> load the protein domains.

############################################################
## Proteins, EnsDb, missing
##
## o Fetch all proteins from the EnsDb database.
## TODO @jo: don't know yet how to handle Uniprot IDs and uniprot table data.
##'
##' @param file EnsDb database from which protein data should/can be retrieved.
##'
##' @param uniprotIds missing.
##'
##' @param loadProteinDomains Logical of length 1 defining whether protein
##' domains within the proteins' sequences should be loaded too. If \code{TRUE},
##' protein domains are loaded and added into the \code{@pranges} slot of the
##' object.
##'
##' @param filter Single object extending
##' \code{\linkS4Object[ensembldb]{BasicFilter}} or \code{list} of such objects
##' to retrieve only specific proteins from the database. See
##' \code{\link[ensembldb]{proteins}} for more information.
##'
##' @param columns Character vector specifying the columns that should be
##' returned from the database. By default, all columns from the \emph{protein}
##' database table are returned. Use the \code{\link[ensembldb]{listColumns}}
##' function to get a list of all supported columns. Note that exon-related
##' columns are not supported for the \code{Proteins} method.
##'
##' @param fetchLRG Logical indicating whether proteins for Locus Reference
##' Genes (LRG) should be retrieved too. By default LRG genes will not be
##' fetched as they do have a 1:n mapping between transcripts and proteins.
##'
##' @param ... Additional arguments to be passed to the
##' \code{\link[ensembldb]{proteins}} method that is used to fetch the data.
##'
##' @noRd
setMethod("Proteins",
          signature(file = "EnsDb", uniprotIds = "missing"),
          function(file, uniprotIds, loadProteinDomains = TRUE,
                   filter = list(), columns = NULL, fetchLRG = FALSE, ...) {
              ## Get the data from EnsDb.
              if (!hasProteinData(file))
                  stop("The provided 'EnsDb' does not contain protein annotations!")
              ## Check 'columns':
              ## o add all columns from the protein table.
              ## o ensure we have no column from the exon or tx2exon table.
              ## o load, if necessary, columns from the protein_domain table.
              columns <- unique(c(columns, listColumns(file, "protein")))
              not_allowed <- listColumns(file, c("exon", "tx2exon"))
              not_allowed <- not_allowed[not_allowed != "tx_id"]
              if (any(columns %in% not_allowed)) {
                  warning("For 'Proteins', fetching exon-related columns is not",
                          " allowed! Columns ",
                          paste0("'", columns[columns %in% not_allowed], "'",
                                 collapse = ", "), " have been removed.")
                  columns <- columns[!(columns %in% not_allowed)]
              }
              if (loadProteinDomains) {
                  columns <- unique(c(columns,
                                      listColumns(file, "protein_domain")))
              }
              if (!fetchLRG) {
                  ##
              }
              ## Now fetch the data:
              res <- ensembldb::proteins(file, filter = filter,
                                         columns = columns,
                                         return.type = "data.frame")
              if (nrow(res) == 0)
                  return(Proteins())
              ## Get the unique data for protein:
              ## For now we're removing the Uniprot IDs. Problem is the 1:n
              ## mapping from protein_id to uniprot_id. Solutions:
              ## 1) return a redundant list of proteins, i.e. replicate the
              ##    proteins, one for each Uniprot ID.
              ## 2) return a unique list of proteins adding the Uniprot IDs as
              ##    a list to the mcols of the AAStringSet (i.e. acols).
              pids <- unique(res$protein_id)
              prt_cn <- colnames(res)[!(colnames(res) %in%
                                        listColumns(file, c("protein_domain",
                                                            "uniprot")))]
              prt <- unique(res[, c("protein_id", prt_cn)])
              aa <- AAStringSet(prt$protein_sequence)
              names(aa) <- prt$protein_id
              mcols(aa) <- prt[, !colnames(prt) %in% c("protein_id",
                                                       "protein_sequence"),
                               drop = FALSE]
              ## Fetch protein domain data
              if (loadProteinDomains) {
                  prng <- IRangesList(list())
                  prt_dom_cols <- listColumns(file, "protein_domain")
                  prt_dom_mcols <- prt_dom_cols[!(prt_dom_cols %in%
                                                  c("protein_domain_id",
                                                    "prot_dom_start",
                                                    "prot_dom_end"))]
                  prt_dom <- unique(res[, prt_dom_cols, drop = FALSE])
                  ## Remove empty ones
                  prt_dom <- prt_dom[!is.na(prt_dom$prot_dom_start), ,
                                     drop = FALSE]
                  if (nrow(prt_dom) > 0) {
                      prng <- IRanges(start = prt_dom$prot_dom_start,
                                      end = prt_dom$prot_dom_end)
                      names(prng) <- prt_dom$protein_domain_id
                      mcols(prng) <- prt_dom[, prt_dom_mcols]
                      prng <- split(prng,
                                    f = factor(prt_dom$protein_id,
                                               levels = unique(prt_dom$protein_id)))
                  }
                  ## Fill empty ones.
                  if (!all(pids %in% prt_dom$protein_id)) {
                      miss_prt <- pids[!(pids %in% prt_dom$protein_id)]
                      tmp <- replicate(length(miss_prt),
                                       .emptyProtDomIRanges(cols = prt_dom_mcols))
                      names(tmp) <- miss_prt
                      prng <- c(prng, IRangesList(tmp))
                  }
                  ## Ensure correct ordering.
                  prng <- prng[pids]
              } else {
                  prng <- IRangesList(replicate(length(pids), IRanges()))
                  names(prng) <- pids
              }
              ## Create the Proteins object.
              pr <- new("Proteins", aa = aa, pranges = prng,
                        metadata = list(created = date(),
                                        ensembl_version = ensemblVersion(file)))
              return(pr)
          })
## Simple helper function to create an IRanges with mcols.
.emptyProtDomIRanges <- function(cols) {
    res <- IRanges()
    tmp <- DataFrame(matrix(ncol = length(cols), nrow = 0))
    colnames(tmp) <- cols
    mcols(res) <- tmp
    return(res)
}


setMethod("[", "Proteins",
          function(x, i, j = "missing", ..., drop) {
              if (!missing(j) || length(list(...)) > 0L)
                  stop("invalid subsetting")
              if (missing(i) || (is.logical(i) && all(i)))
                  return(x)
              if (is.logical(i))
                  i <- which(i)
              if (is.character(i))
                  i <- match(i, seqnames(x))
              if (!is.numeric(i) || any(is.na(i)))
                  stop("invalid subsetting")
              if (any(i < 1) || any(i > length(x)))
                  stop("subscript out of bounds")

              x@aa <- x@aa[i]
              x@pranges <- x@pranges[i]
              x
          })


## accessors
setMethod("pfeatures", "Proteins",
          function(x) extractAt(aa(x), unname(pranges(x))))

setMethod("pranges", "Proteins",
          function(x) x@pranges)

setReplaceMethod("pranges",
                 c("Proteins", "CompressedIRangesList"),
                 function(object, value)
                     replacePranges(object, value))

setReplaceMethod("acols",
                 c("Proteins", "DataFrame"),
                 function(object, value)
                     replaceAcols(object, value))

setMethod("length", "Proteins",
          function(x) length(x@aa))

setMethod("metadata", "Proteins",
          function(x) x@metadata)

## signature is (x, ..., value)
setReplaceMethod("metadata", "Proteins",
                 function(x, name, value) {
                     if (name == "created")
                         stop("Creation date can't be modified.")
                   x@metadata[[name]] <- value
                   x
               })

setMethod("pmetadata", "Proteins",
          function(x) {
              if (!is.null(x@pranges@unlistData@elementMetadata)) {
                  SplitDataFrameList(lapply(x@pranges, mcols))
              } else {
                  NULL
              }
          })

setMethod("ametadata", "Proteins",
          function(x) mcols(x@aa))

setMethod("seqnames","Proteins",
          function(x) names(aa(x)))


setMethod("avarLabels", "Proteins",
          function(object) names(aa(object)@elementMetadata))

setMethod("pvarLabels", "Proteins",
          function(object) {
              if (!isEmpty(pranges(object))) {
                names(pranges(object)@unlistData@elementMetadata@listData)
              } else {
                NULL
              }
          })

setMethod("[[", "Proteins",
          function(x, i, j = missing, ..., drop = TRUE) x@aa[[i]])

setMethod("aa", "Proteins", function(x) x@aa)


## Methods
setMethod("addIdentificationData",
          c("Proteins", "character"),
          function(object, id, rmEmptyRanges = TRUE, par = Pparams()) {
              .addIdentificationDataProteins(object, filenames = id,
                                             rmEmptyRanges = rmEmptyRanges,
                                             par = par)
          })

setMethod("addPeptideFragments",
          c("Proteins", "character"),
          function(object, filenames, rmEmptyRanges = TRUE, par = Pparams()) {
              .addPeptideFragmentsProteins(object, filenames = filenames,
                                           rmEmptyRanges = rmEmptyRanges,
                                           par = par)
          })

setMethod("cleave", "Proteins",
          function(x, enzym = "trypsin", missedCleavages = 0, ...) {
              rng <- cleavageRanges(x = x@aa, enzym = enzym,
                                          missedCleavages = missedCleavages,
                                          ...)
              x@pranges <- IRangesList(lapply(rng, function(r) {
                  mc <- missedCleavages[cumsum(start(r) == 1L)]
                  mcols(r) <- DataFrame(MissedCleavages = Rle(mc))
                  r
              }))
              x
          })

setMethod("plot",
          signature(x = "Proteins", y = "missing"),
          function(x, y, ...) .plotProteins(x, ...))

setMethod("pfilter",
          "Proteins",
          function(x, mass = NULL, len = NULL, ...) {
              .pfilterProteins(x, mass = mass, len = len, ...)
          })

setMethod("show", "Proteins",
          function(object) {
              topics <- c("S4 class type",
                          "Class version",
                          "Created",
                          "Number of Proteins")
              topics <- format(topics, justify = "left")
              n <- length(object)
              values <- c(class(object),
                          tail(as(classVersion(object), "character"), 1L),
                          object@metadata$created,
                          n)

              cat(paste0(topics, ": ",  values, collapse = "\n"), sep = "\n")
              if (length(object) > 0) {
                  sn <- seqnames(object)
                  ln <- length(object)
                  cat("Sequences:\n  "); htcat(sn, n = 2)
                  cat("Sequence features:\n  "); htcat(avarLabels(object), n = 2)
                  cat("Peptide features:")
                  if (isEmpty(pranges(object))) cat(" None\n")
                  else {
                      cat("\n  ")
                      htcat(pvarLabels(object), n = 2)
                  }
              }
          })


## New suggested methods:
## + width: return the width(x@aa)
## + names,Proteins: returns seqnames,Proteins

## internal use only; not exported

setMethod("aaranges",
          "Proteins",
          function(x, ...) {
              .aarangesProteins(x, ...)
          })
