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
    "Chromosome", "..maf_cols", "..maf_cols2",
    "Sample", "Tumor_Sample_Barcode", ":=",

    # seg_to_arm_cnv
    "width", "arm_width", "seqnames", "arm",

    # stat_between_test_helper
    "term", "method", "term_x", "term_y",

    # stat_cor_test
    "value",

    # tcga_get_cli_biotab, tcga_get_cli_indexed, tcga_get_cli_xml,
    # tcga_remove_duplicated_samples
    "bcr_patient_uuid", "bcr_sample_barcode", "sample_type",
    "bcr_patient_barcode","sample_barcode", "sample_vial_barcode",
    "submitter_id", "analyte", "plate", "portion"

  )
)
