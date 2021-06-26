.onLoad <- function(libname, pkgname) {

  version <- utils::packageDescription(pkgname, fields = "Version")

  packageStartupMessage(
    pkgname, " version ", version, "\n",
    "######################################\n",
    "## ##### ####### ### ######    #######\n",
    "### ### ###### ### ### ####### #######\n",
    "##### ######## ####### ####### #######\n",
    "##### ######### ##### ##### ## #######\n",
    "##### ###########   ########  ########\n",
    "######################################\n",
    "Miscellaneous fns for Yun and Jin!!"
  )
  invisible()

}

.onAttach <- function(libname, pkgname) {

}
