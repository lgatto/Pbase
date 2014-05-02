.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nThis is Pbase version ",
    utils::packageVersion("Pbase"), "\n")
}
