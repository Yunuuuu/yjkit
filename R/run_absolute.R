#' Automate ABSOLUTE calling for multiple samples in parallel ways
#'
#' @description This function is based on package
#'   \href{https://github.com/ShixiangWang/DoAbsolute}{DoAbsolute} by adding
#'   parallel evaluation and adjusting personal convention
#'
#' @details \href{https://www.nature.com/articles/nbt.2203}{ABSOLUTE} is a
#'   famous software developed by Broad Institute, However, the
#'   \code{\link[ABSOLUTE]{RunAbsolute}} function points to estimate one sample
#'   each time and sets no default values. \code{\link{run_absolute}} helps
#'   users set default parameters based on
#'   \href{https://www.genepattern.org/modules/docs/ABSOLUTE}{ABSOLUTE
#'   documentation} and provides an uniform interface to input data easily and
#'   run \code{\link[ABSOLUTE]{RunAbsolute}} parallelly.
#'
#'   If calling for a sample failed, the error message will be written to
#'   \code{error.log} under subdirectory \code{RunAbsolute} of
#'   \code{results_dir} directory.
#'
#'   More detail about how to analyze ABSOLUTE results please see
#'   \href{https://www.genepattern.org/analyzing-absolute-data}{analyzing-absolute-data}.
#'
#' @param seg a \code{data.frame} containing columns "Sample", "Chromosome",
#'   "Start", "End", "Num_Probes", "Segment_Mean".
#' @param maf MAF, default is \code{NULL}, can provided as \code{data.frame}.
#' @param sigma_p Provisional value of excess sample level variance used for
#'   mode search. Default: \code{0}
#' @param max_sigma_h Maximum value of excess sample level variance (Eq. 6).
#'   Default: \code{0.015}
#' @param min_ploidy Minimum ploidy value to consider. Solutions implying lower
#'   ploidy values will be discarded. Default: \code{0.95}
#' @param max_ploidy Maximum ploidy value to consider. Solutions implying
#'   greater ploidy values will be discarded. Default: \code{10}
#' @param primary_disease Primary disease of the sample. Default: \code{NA}
#' @param platform one of \code{"SNP_6.0"}, \code{"Illumina_WES"},
#'   \code{"SNP_250K_STY"}. Default: \code{"SNP_6.0"}
#' @param results_dir directory path used to store result files. Default:
#'   \code{here::here("results", "ABSOLUTE")}
#' @param max_as_seg_count Maximum number of allelic segments. Samples with a
#'   higher segment count will be flagged as 'failed'. Default: \code{1500}
#' @param max_non_clonal Maximum genome fraction that may be modeled as
#'   non-clonal (subclonal SCNA). Solutions implying greater values will be
#'   discarded. Default: \code{0.05}
#' @param max_neg_genome Maximum genome fraction that may be modeled as
#'   non-clonal with copy-ratio below that of clonal homozygous deletion.
#'   Solutions implying greater values will be discarded. Default: \code{0.005}
#' @param copy_num_type The type of copy number to be handled. Either total or
#'   allelic. total is what this code for. Default: \code{"total"}
#' @param min_mut_af Minimum mutation allelic fraction. Mutations with lower
#'   allelic fractions will be filtered out before analysis. Default: \code{0.1}
#' @param verbose if \code{TRUE}, print extra info. Default: \code{FALSE}
#' @param BPPARAM Default: \code{\link[BiocParallel:bpparam]{bpparam()}}. \cr
#'   \code{NULL} means \code{switch (Sys.info()[["sysname"]],} \cr
#'   \code{Darwin = ,} \cr \code{Linux = BiocParallel::MulticoreParam(),} \cr
#'   \code{Windows = BiocParallel::SnowParam()) }
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @return  Side effect. All ABSOLUTE called results (see
#'   \code{\link[ABSOLUTE]{RunAbsolute}}) were kept in directory
#'   \code{file.path(results_dir, "RunAbsolute")}, all summarized multiple
#'   ABSOLUTE results (see \code{\link[ABSOLUTE]{CreateReviewObject}}) were kept
#'   in \code{file.path(results_dir, "CreateReviewObject")}, and all reviewed
#'   results were kept in \code{ file.path(results_dir, "reviewed")}
#' @examples
#' \donttest{
#' seg <- readRDS(
#' system.file("extdata", "run_absolute_example_seg.rds", package = "yjtools")
#' )
#'
#' maf <- readRDS(
#' system.file("extdata", "run_absolute_example_maf.rds", package = "yjtools")
#' )
#'
#' run_absolute(seg = seg, maf = maf,
#' results_dir = file.path(tempdir(),"results", "ABSOLUTE"),
#' BPPARAM = NULL
#' )
#' }
#' @references Carter, S., Cibulskis, K., Helman, E. et al. Absolute
#'   quantification of somatic DNA alterations in human cancer. Nat Biotechnol
#'   30, 413–421 (2012). \url{https://doi.org/10.1038/nbt.2203}
#' @import data.table
#' @export
run_absolute <- function(
  seg, maf = NULL, sigma_p = 0, max_sigma_h = 0.015, min_ploidy = 0.95,
  max_ploidy = 10, primary_disease = NA,
  platform = c("SNP_6.0", "Illumina_WES", "SNP_250K_STY"),
  results_dir = here::here("results", "ABSOLUTE"),
  max_as_seg_count = 1500, max_neg_genome = 0.005,
  max_non_clonal = 0.05, copy_num_type = c("total", "allelic"),
  min_mut_af = 0.1,  verbose = FALSE,
  BPPARAM = BiocParallel::bpparam()){

  if (!requireNamespace("ABSOLUTE", quietly = TRUE)){
    stop("ABSOLUTE needed for this function to work. Please install it",
         call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)){
    stop("data.table needed for this function to work. Please install it",
         call. = FALSE)
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE)){
    stop("BiocParallel needed for this function to work. Please install it",
         call. = FALSE)
  }

  if (verbose) message("Setting results directory as",
                       results_dir, appendLF = TRUE)
  dir_create2(results_dir)

  temp_dir <- file.path(tempdir(), "ABSOLUTE")
  if (verbose) message("Setting temporary directory as ",
                       temp_dir, appendLF = TRUE)
  dir_create2(temp_dir)
  on.exit( unlink(temp_dir, recursive = TRUE, force = TRUE),
           add = TRUE)

  if (verbose) message("Confirming seg (and maf) data with right patterns.",
                       appendLF = TRUE)
  absolute_data <- run_absolute_validate_seg_and_maf_data(seg = seg, maf = maf)

  # create temp files to generate input for RunAbsolute

  if (verbose) message("Writing each seg (and maf) data of individual into ",
                       "different files based on samples column in seg data.",
                       appendLF = TRUE)

  absolute_filepath <- run_absolute_prepare_seg_and_maf_data(
    seg = absolute_data[["seg"]],
    maf = absolute_data[["maf"]],
    temp_dir = temp_dir
  )

  # match options
  platform <- match.arg(platform)
  copy_num_type <- match.arg(copy_num_type)

  if (length(absolute_filepath[["sample_id"]]) > 0) {

    if (verbose) message("Running ABSOLUTE algorithm...", appendLF = TRUE)

    if (is.null(BPPARAM)) BPPARAM <- switch (
      Sys.info()[["sysname"]],
      Darwin = ,
      Linux = BiocParallel::MulticoreParam(progressbar = verbose),
      Windows = BiocParallel::SnowParam(progressbar = verbose)
    )

    run_absolute_dir <- file.path(
      results_dir, "RunAbsolute"
    )

    BiocParallel::bplapply(
      absolute_filepath[["sample_id"]],
      function(i, seg_filepath, maf_filepath, results_dir,
               sigma_p, max_sigma_h, min_ploidy, max_ploidy,
               primary_disease, platform, max_as_seg_count,
               max_non_clonal, max_neg_genome, copy_num_type, min_mut_af,
               verbose){

        suppressPackageStartupMessages( loadNamespace("ABSOLUTE") )

        if (is.na(maf_filepath[[i]])) {
          maf_fn <- NULL
        } else {
          maf_fn <- maf_filepath[[i]]
        }

        safe_runAbsolute(
          seg_dat_fn = seg_filepath[i],
          maf_fn = maf_fn, sample_name = i,
          sigma_p = sigma_p, max_sigma_h = max_sigma_h,
          min_ploidy = min_ploidy, max_ploidy = max_ploidy,
          primary_disease = primary_disease, platform = platform,
          results_dir = results_dir, max_as_seg_count = max_as_seg_count,
          max_non_clonal = max_non_clonal, max_neg_genome = max_neg_genome,
          copy_num_type = copy_num_type, min_mut_af = min_mut_af,
          verbose = verbose
        )

      },
      seg_filepath = absolute_filepath[["seg"]],
      maf_filepath = absolute_filepath[["maf"]],
      sigma_p = sigma_p, max_sigma_h = max_sigma_h,
      min_ploidy = min_ploidy, max_ploidy = max_ploidy,
      primary_disease = primary_disease, platform = platform,
      results_dir = run_absolute_dir, max_as_seg_count = max_as_seg_count,
      max_non_clonal = max_non_clonal, max_neg_genome = max_neg_genome,
      copy_num_type = copy_num_type, min_mut_af = min_mut_af,
      verbose = verbose, BPPARAM = BPPARAM)

    if (verbose) message("RunAbsolute done.", appendLF = TRUE)

    run_absolute_files <- file.path(
      run_absolute_dir,
      paste0(absolute_filepath[["sample_id"]], ".ABSOLUTE.RData")
    )

    run_absolute_files <- run_absolute_files[
      file.exists(run_absolute_files)
    ]

    if (length(run_absolute_files) == 0) {
      stop("No RunAbsolute results file to proceed.",
           call. = FALSE)
    }

    summarize_dir <- file.path(results_dir, "CreateReviewObject")

    if (dir.exists(summarize_dir)) {
      if (verbose) message("Removing previous summarize results directory.",
                           appendLF = TRUE)
      unlink(summarize_dir, recursive = TRUE)
    }

    if (verbose) message("Summarizing multiple ABSOLUTE runs...",
                         appendLF = TRUE)

    suppressWarnings(ABSOLUTE::CreateReviewObject(
      obj.name = "SummarizeAbsolute",
      absolute.files = run_absolute_files,
      indv.results.dir = summarize_dir,
      copy_num_type = copy_num_type,
      plot.modes = TRUE,
      verbose = verbose
    ))

    if (verbose) message("Absolute summarize done.\nPrepare auto-reviewing...",
                         appendLF = TRUE)

    pp_call_fn <- file.path(summarize_dir,
                            "SummarizeAbsolute.PP-calls_tab.txt")
    modes_fn <- file.path(summarize_dir,
                          "SummarizeAbsolute.PP-modes.data.RData")

    suppressWarnings(ABSOLUTE::ExtractReviewedResults(
      reviewed.pp.calls.fn = pp_call_fn,
      analyst.id = "YJ",
      modes.fn = modes_fn,
      out.dir.base = results_dir,
      obj.name = "ReviewAbsolute",
      copy_num_type = copy_num_type,
      verbose = verbose
    ))

    if (verbose) message("Reviewing ABSOLUTE summary done.",
                         appendLF = TRUE)

  }

  cat("ABSOLUTE algorithm Done.")
}


# run_absolute utility functions ----------------------------------------------

safe_runAbsolute <- function(
  seg_dat_fn, maf_fn, sample_name, sigma_p, max_sigma_h,
  min_ploidy, max_ploidy, primary_disease, platform,
  results_dir, max_as_seg_count, max_non_clonal,
  max_neg_genome, copy_num_type, min_mut_af,
  verbose){

  error_log_file <- file.path(results_dir, "error.log")

  if (file.exists(error_log_file)) {

    message("Removing previous error log file in ", results_dir,
            appendLF = TRUE)
    if ( !file.remove(error_log_file) ){
      warning("Cannot remove previous error log file.", call. = FALSE)
    }

  }

  tryCatch(
    {
      suppressWarnings(ABSOLUTE::RunAbsolute(
        seg.dat.fn = seg_dat_fn, maf.fn = maf_fn,
        sample.name = sample_name,
        sigma.p = sigma_p, max.sigma.h = max_sigma_h,
        min.ploidy = min_ploidy, max.ploidy = max_ploidy,
        primary.disease = primary_disease, platform = platform,
        results.dir = results_dir, max.as.seg.count = max_as_seg_count,
        max.non.clonal = max_non_clonal,
        max.neg.genome = max_neg_genome,
        copy_num_type = copy_num_type,
        min.mut.af = min_mut_af,
        verbose = verbose
      ))
    },
    error = function(cnd) {

      sink(error_log_file, append = TRUE)
      on.exit(sink(file = NULL), add = TRUE)

      cat("Detecting error in sample: ", sample_name, "\n")

      cat("Error message:\n", conditionMessage(cnd), "\n")

      if (grepl("mutations left", conditionMessage(cnd))) {
        cat("Try to fix error by removing maf file.\n")
        tryCatch(
          {
            suppressWarnings(ABSOLUTE::RunAbsolute(
              seg.dat.fn = seg_dat_fn,  maf.fn = NULL,
              sample.name = sample_name,
              sigma.p = sigma_p, max.sigma.h = max_sigma_h,
              min.ploidy = min_ploidy, max.ploidy = max_ploidy,
              primary.disease = primary_disease, platform = platform,
              results.dir = results_dir,
              max.as.seg.count = max_as_seg_count,
              max.non.clonal = max_non_clonal,
              max.neg.genome = max_neg_genome,
              copy_num_type = copy_num_type,
              verbose = verbose
            ))
            cat("Fixing successfully!\n")
          },
          error = function(cnd) {
            cat("Fixing failed. Skipping this sample.\n")
          }
        )
      } else {
        cat("Skipping this sample.\n")
      }

      cat("========\n\n")
    }
  )

}

dir_create2 <- function(x) {

  if (!dir.exists(x)) {
    message("Directory ", x, " does not exists. try to create it.")

    if(!dir.create(x, recursive = TRUE)){
      stop("Cannot create directory ", x,
           call. = FALSE)
    }

  }

}

# validate seg and maf data to have corresponding columns -----------------

run_absolute_validate_seg_and_maf_data <- function(seg, maf = NULL){

  if (!inherits(seg, "data.frame")){
    stop("The class of seg must inherited from data.frame including ",
         "data.frame, data.table, or tibble, but not ",
         class(seg), call. = FALSE)
  }

  seg <- data.table::as.data.table(seg)

  # check seg data
  seg_cols <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

  if (!all(seg_cols %in% colnames(seg))) {

    stop("Object seg should include columns ",
         paste(seg_cols[-length(seg_cols)], collapse = ", "), " ",
         "and ", seg_cols[length(seg_cols)], "; ",
         "we couldn't find columns ",
         paste(setdiff(seg_cols, colnames(seg)), collapse = ", "), ".",
         call. = FALSE)

  }

  seg <- seg[, seg_cols, with = FALSE]

  seg[, Chromosome := as.character(Chromosome)]
  seg[, Chromosome := gsub(pattern = "chr", replacement = "",
                           Chromosome, ignore.case = TRUE)]
  seg[, Chromosome := gsub(pattern = "X", replacement = "23",
                           Chromosome, ignore.case = TRUE)]

  seg <- seg[Chromosome %in% as.character(1:23), ]

  # check maf data
  if (!is.null(maf)) {

    if (!inherits(maf, "data.frame")){
      stop("The class of maf must inherited from data.frame including ",
           "data.frame, data.table, or tibble, but not ",
           class(maf), call. = FALSE)
    }
    maf <- data.table::as.data.table(maf)

    maf_cols <- c(
      Tumor_Sample_Barcode = "Tumor_Sample_Barcode",
      Chromosome = "Chromosome",
      Hugo_Symbol = "Hugo_Symbol",
      dbSNP_Val_Status = "dbSNP_Val_Status"
    )

    lack_cols <- setdiff(
      names(maf_cols), colnames(maf)
    )

    if ("Start_position" %in% colnames(maf)) {
      maf_cols <- c(maf_cols, Start_position = "Start_position")
    } else if ("Start_Position" %in% colnames(maf)) {
      maf_cols <- c(maf_cols, Start_Position = "Start_position")
    } else {

      lack_cols <- c(lack_cols, "Start_position or Start_Position")

    }

    if ("t_ref_count" %in% colnames(maf)) {

      maf_cols <- c(maf_cols, t_ref_count = "t_ref_count")

    } else if ("i_t_ref_count" %in% colnames(maf)) {

      maf_cols <- c(maf_cols, i_t_ref_count = "t_ref_count")

    } else {

      lack_cols <- c(lack_cols, "t_ref_count or i_t_ref_count")

    }

    if ("t_alt_count" %in% colnames(maf)) {
      maf_cols <- c(maf_cols, t_alt_count = "t_alt_count")
    } else if ("i_t_alt_count" %in% colnames(maf)) {
      maf_cols <- c(maf_cols, i_t_alt_count = "t_alt_count")
    } else {

      lack_cols <- c(lack_cols, "t_alt_count or i_t_alt_count")

    }

    if (length(lack_cols) > 0){
      stop("Object maf should include columns ",
           "Tumor_Sample_Barcode, Chromosome, Hugo_Symbol, dbSNP_Val_Status ",
           "Start_position or Start_Position ",
           "t_ref_count or i_t_ref_count ",
           "and t_alt_count or i_t_alt_count. ",
           "we couldn't find columns ",
           paste(lack_cols, collapse = ", "), ".",
           call. = FALSE)
    }

    maf <- maf[, names(maf_cols), with = FALSE]
    colnames(maf) <- maf_cols

    maf[, Chromosome := as.character(Chromosome)]
    maf[, Chromosome := gsub(pattern = "chr", replacement = "", Chromosome, ignore.case = TRUE)]
    maf[, Chromosome := gsub(pattern = "X", replacement = "23", Chromosome, ignore.case = TRUE)]
    maf <- maf[Chromosome %in% as.character(1:23), ]

  }

  list(
    seg = seg,
    maf = maf
  )

}

run_absolute_prepare_seg_and_maf_data <- function(seg, maf = NULL, temp_dir){

  .sample_id. <- unique( seg$Sample )

  if (any(is.na(.sample_id.))){
    warning("NA value found in Sample id in seg data, ",
            "we'll remove it",
            call. = FALSE)
    .sample_id. <- .sample_id.[!is.na(.sample_id.)]

    if (length(.sample_id.) == 0) stop(
      "We couldn't find any samples (excluding NA values)"
    )
  }

  # prepare seg data
  seg <- seg[Sample %in% .sample_id., ]
  .seg_filepath. <- file.path(temp_dir,
                              paste0(.sample_id., ".seg"))

  names(.seg_filepath.) <- .sample_id.

  seg[, data.table::fwrite(
    x = .SD,
    file = .seg_filepath.[[unlist(.BY)]],
    sep = "\t"
  ), by = Sample]


  # prepare maf data

  if (!is.null(maf)){

    maf <- maf[Tumor_Sample_Barcode %in% .sample_id., ]
    maf[, .group_col. := Tumor_Sample_Barcode]
    .maf_filepath. <- ifelse(
      .sample_id. %in% maf[["Tumor_Sample_Barcode"]],
      file.path(temp_dir, paste0(.sample_id., ".maf")),
      NA_character_
    )

    names(.maf_filepath.) <- .sample_id.

    maf[, data.table::fwrite(
      x = .SD,
      file = .maf_filepath.[[unlist(.BY)]],
      sep = "\t"
    ), by = .group_col.]

  }

  list(
    sample_id = .sample_id.,
    seg = .seg_filepath.,
    maf = .maf_filepath.
  )
}
