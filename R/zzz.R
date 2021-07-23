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
  c(# ggplot2 utils functions
    "x", "y", "paletteer",

    # stats
    "stats",

    # run_absolute
    "Chromosome", "Sample", "Tumor_Sample_Barcode",
    ":=", ".group_col.",

    # seg_to_arm_cnv
    "width", "arm_width", "seqnames", "arm", "CNV", "seg_length_frac",

    # stat_between_test_helper
    "term", "method", "term_x", "term_y",

    # stat_cor_test
    "value"

  )
)
