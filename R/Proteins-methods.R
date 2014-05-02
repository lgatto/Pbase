setMethod("Proteins",
          signature(file = "missing", uniprotIds = "character"),
          function(file, uniprotIds, ...) {
            .toBeImplemented()
          })

setMethod("[", "Proteins",
          function(x, i, j, ..., drop = FALSE) {
            x@aa <- x@aa[i, ]
            return(x)
          })

setMethod("[[", "Proteins",
          function(x, i, j, ..., drop = TRUE) {
            return(x@aa[[i]])
          })

setMethod("aa",
          "Proteins",
          function(x)x@aa)

setMethod("ametadata",
          "Proteins",
          function(x)mcols(x@aa))

setMethod("cleave", "Proteins",
          function(x, ...) {
              x@pfeatures <- cleave(aa(x), ...)
              x
          })

setMethod("length",
          "Proteins",
          function(x)length(x@aa))

setMethod("metadata",
          "Proteins",
          function(x)x@metadata)

setMethod("plot",
          "Proteins",
          function(x, ...).plotProteins(x, ...))

setMethod("pmetadata",
          "Proteins",
          function(x)mcols(x@pfeatures))

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

  cat("Sequences:", tail(capture.output(object@aa), -1), sep = "\n")
})

