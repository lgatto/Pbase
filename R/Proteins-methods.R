setMethod("Proteins",
          signature(file = "character", uniprotIds = "missing"),
          function(file, uniprotIds, ...) {
            .ProteinsFromFasta(filenames = file, ...)
          })

setMethod("Proteins",
          signature(file = "missing", uniprotIds = "character"),
          function(file, uniprotIds, ...) {
            .toBeImplemented()
          })

setMethod("[", "Proteins",
          function(x, i, j, ..., drop = FALSE) {
            x@mcols <- x@mcols[i, , drop = drop ]
            x@seq <- x@seq[i, ]
            return(x)
          })

setMethod("[[", "Proteins",
          function(x, i, j, ..., drop = TRUE) {
            return(x@seq[[i]])
          })

setMethod("length",
          "Proteins",
          function(x)length(x@seq))

setMethod("mcols",
          "Proteins",
          function(x)x@mcols)

setMethod("metadata",
          "Proteins",
          function(x)x@metadata)

setMethod("show",
          "Proteins",
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

  cat("Sequences:", tail(capture.output(object@seq), -1), sep = "\n")
})

