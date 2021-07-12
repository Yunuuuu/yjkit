## code to prepare `DATASET` dataset goes here

# run_arm_cnv --------------------------------------------------------

# ** ref_cytoband ---------------------------------------------------------

anno_hub <- AnnotationHub::AnnotationHub(localHub = FALSE)
AnnotationHub::query(
  anno_hub, c("UCSC", "Homo sapiens", "cytoband"),
  ignore.case=TRUE
)
# AH53177 | UCSC cytoBand track for hg19
# AH53178 | UCSC cytoBand track for hg38

run_arm_cnv_ref_cytoband_hg19 <- anno_hub[["AH53177"]]
run_arm_cnv_ref_cytoband_hg38 <- anno_hub[["AH53178"]]

# run_cibersort -----------------------------------------------------------

run_cibersort_lm22 <- readr::read_tsv(
  here::here("data-raw", "LM22.txt"),
  col_names = TRUE
)

usethis::use_data(
  run_arm_cnv_ref_cytoband_hg19,
  run_arm_cnv_ref_cytoband_hg38,
  run_cibersort_lm22,
  internal = TRUE,
  overwrite = TRUE,
  compress = "gzip",
  version = 3
)
