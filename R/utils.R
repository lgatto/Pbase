.toBeImplemented <- function() {
  parentCall <- sys.call(-1L)
  stop(deparse(parentCall), " is not implemented yet!", call.=FALSE)
}
