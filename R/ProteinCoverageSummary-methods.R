setMethod("show",
          "ProteinCoverageSummary",
          function(object) {

  callNextMethod(object)

  cv <- .calculateProteinCoverageSummary(object)
  overall <- attr(cv, "overall")
  nm <- names(cv)

  cv <- sprintf("%0.3f", cv)

  n <- length(object)

  if (n > 10L) {
    nm <- c(nm[1L:5L], "...", nm[6L:10L])
    cv <- c(cv[1L:5L], "...", cv[6L:10L])
  }

  nm <- format(nm, justify = "left")

  cat("Summary (overall coverage ", sprintf("%0.3f", overall), "):\n", sep="")
  cat(paste(nm, cv, collapse = "\n"), sep = "\n")
})

