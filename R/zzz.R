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

if(getRversion() >= "2.15.1") utils::globalVariables(
  c(# run_absolute
    "Chromosome", "Sample", "Tumor_Sample_Barcode",
    ":=", ".group_col."
  )
)
