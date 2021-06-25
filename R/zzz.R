.onLoad <- function(libname, pkgname) {

  version <- utils::packageDescription(pkgname, fields = "Version")

  # packageStartupMessage(
  #
  # )
  invisible()

}

.onAttach <- function(libname, pkgname) {

}
