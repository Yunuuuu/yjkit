.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  packageStartupMessage(
    "Package: ", pkgname, " (version: ", version, ")\n",
    "######################################\n",
    "##### ##### ###### ### ######    #####\n",
    "###### ### ##### ### ### ####### #####\n",
    "######## ####### ####### ####### #####\n",
    "######## ######## ##### ##### ## #####\n",
    "######## ##########   ########  ######\n",
    "######################################\n",
    "Miscellaneous fns for Yun and Jin!!"
  )
  invisible()
}
