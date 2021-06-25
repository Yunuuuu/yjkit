.onLoad <- function(libname, pkgname) {

  version <- packageDescription(pkgname, fields = "Version")

  # packageStartupMessage(
  #
  # )
  invisible()

}

.onAttach <- function(libname, pkgname) {

}
