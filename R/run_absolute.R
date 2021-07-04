###################################################
######## Automate ABSOLUTE Calling ################
###################################################
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
#'   \strong{error.log} under the \code{results_dir} directory.
#'
#'   More detail about how to analyze ABSOLUTE results please see
#'   \href{https://www.genepattern.org/analyzing-absolute-data}{analyzing-absolute-data}.
#'
#' @param seg a \code{data.frame} or a \code{file path} read as a data.fram
#'   containing columns "Sample", "Chromosome", "Start", "End", "Num_Probes",
#'   "Segment_Mean".
#' @param maf MAF, default is \code{NULL}, can provided as \code{data.frame} or
#'   \code{file path}.
#' @param sigma_p Provisional value of excess sample level variance used for
#'   mode search. Default: \code{0}
#' @param max_sigma_h Maximum value of excess sample level variance (Eq. 6).
#'   Default: \code{0.015}
#' @param min_ploidy Minimum ploidy value to consider. Solutions implying lower
#'   ploidy values will be discarded. Default: \code{0.95}
#' @param max_ploidy Maximum ploidy value to consider. Solutions implying
#'   greater ploidy values will be discarded. Default: \code{10}
#' @param primary_disease Primary disease of the sample. Default: \code{NA}
#' @param platform one of "SNP_6.0", "Illumina_WES", "SNP_250K_STY". Default:
#'   \code{"SNP_6.0"}
#' @param temp_dir directory path used to store tempory files. Default: ABSOLUTE
#'   subdirectory under \code{tempdir()}
#' @param clean_temp if \code{TRUE}, auto-clean temp dir at the end. Default:
#'   \code{FALSE}
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
#' allelic. total is what this code for. Default: \code{"total"}
#' @param min_mut_af Minimum mutation allelic fraction. Mutations with lower
#'   allelic fractions will be filtered out before analysis. Default: \code{0.1}
#' @param keep_all_results if \code{TRUE}, clean all results, otherwise clean
#'   result directory and keep most important results. Default: \code{TRUE}
#' @param recover if \code{TRUE}, recover previous unfinished work. This is
#'   helpful when program stop unexpectedly when \code{clean.temp} is
#'   \code{FALSE}. Default: \code{FALSE}
#' @param verbose if \code{TRUE}, print extra info. Default: \code{FALSE}
#' @param BPPARAM see \code{\link[BiocParallel:register]{bpparam}}. Default:
#'   \code{NULL} means \cr \code{switch (Sys.info()[["sysname"]],} \cr
#'   \code{Linux = BiocParallel::MulticoreParam(),} \cr \code{Windows =
#'   BiocParallel::SnowParam()) }
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @return All reviewed results (see
#'   \code{\link[ABSOLUTE]{ExtractReviewedResults}}) were in directory
#'   \code{results_dir}.
#'
#'   if keep_all_results was TRUE, all ABSOLUTE called results (see
#'   \code{\link[ABSOLUTE]{RunAbsolute}}) were kept in directory
#'   \code{file.path(results_dir, "call_res")}, all summarized multiple ABSOLUTE
#'   results (see \code{\link[ABSOLUTE]{CreateReviewObject}}) were kept in
#'   \code{file.path(results_dir, "summary_res")}, and all reviewed samples were
#'   kept in \code{ file.path(results_dir, "sample_summaized")}
#' @examples
#' \donttest{
#' seg <- readRDS(
#' system.file("extdata", "run_absolute_example_seg.rds", package = "yjkit")
#' )
#'
#' maf <- readRDS(
#' system.file("extdata", "run_absolute_example_maf.rds", package = "yjkit")
#' )
#'
#' run_absolute(seg = seg, maf = maf,
#' results_dir = file.path(tempdir(),"results", "ABSOLUTE")
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
  temp_dir = file.path(tempdir(), "ABSOLUTE"), clean_temp = FALSE,
  results_dir = here::here("results", "ABSOLUTE"),
  max_as_seg_count = 1500, max_neg_genome = 0.005,
  max_non_clonal = 0.05, copy_num_type = c("total", "allelic"),
  min_mut_af = 0.1, keep_all_results = TRUE, recover = FALSE,
  verbose = FALSE, BPPARAM = NULL){

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

  cat("-> Setting results directory as", results_dir, "\n")
  dir_create2(results_dir)

  if (file.exists(file.path(results_dir, "error.log"))) {
    unlink(file.path(results_dir, "error.log"), recursive = FALSE)
    cat("-> Removed previous error log file.\n")
  }

  cat("-> Setting temp directory as", temp_dir, "\n")

  if (verbose  && !dir.exists(temp_dir)) cat("-> Creating temp directory ...\n")
  dir_create2(temp_dir)
  on.exit(
    if (clean_temp && dir.exists(temp_dir)) { unlink(temp_dir, recursive = TRUE, force = TRUE) }
  )

  cache_dir <- file.path(temp_dir, "cache")
  dir_create2(cache_dir)

  if (verbose) cat("-> Loading segmentation data...\n")

  if (is.character(seg) && length(seg) == 1) {

    if (file.exists(seg)) {
      seg <- data.table::fread(
        input = seg, data.table = TRUE
      )
    } else {
      stop("file ", seg, " does not exist")
    }

  } else {

    if (inherits(seg, "data.frame")) {
      seg <- data.table::as.data.table(seg)
    } else {
      stop("Unsupport Segmentation Format!")
    }

  }

  if (!is.null(maf)) {

    if (verbose) cat("-> Loading maf data...\n")

    if (is.character(maf)) {

      if (file.exists(maf)) {
        maf <- data.table::fread(input = maf)
      } else {
        stop("file ", maf, " does not exist")
      }

    } else {

      if (inherits(maf, "data.frame")) {
        maf <- data.table::as.data.table(maf)
      } else {
        stop("Unsupport maf Format!")
      }
    }

  }

  if (verbose) cat("-> Checking data format of segmentation file...\n")

  #-- check seg data
  seg_cols <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

  if (!all(seg_cols %in% colnames(seg))) {

    stop("Columns ", paste(seg_cols, collapse = " "),
         " should be in seg")

  } else {

    seg <- seg[, c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")]

  }

  if (verbose) cat("-> Keeping only chr 1-23 for CNV data...\n")

  seg[, Chromosome := as.character(Chromosome)]
  seg[, Chromosome := gsub(pattern = "chr", replacement = "", Chromosome, ignore.case = TRUE)]
  seg[, Chromosome := gsub(pattern = "X", replacement = "23", Chromosome, ignore.case = TRUE)]
  seg <- seg[Chromosome %in% as.character(1:23), ]

  #-- check maf data
  if (!is.null(maf)) {

    maf_cols <- c(
      "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_position", "dbSNP_Val_Status", "t_ref_count", "t_alt_count"
    )
    maf_cols2 <- c(
      "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "dbSNP_Val_Status", "t_ref_count", "t_alt_count"
    )

    if (all(maf_cols %in% colnames(maf))) {

      maf <- maf[, ..maf_cols]

    } else if (all(maf_cols2 %in% colnames(maf))) {

      maf <- maf[, ..maf_cols2]
      colnames(maf) <- maf_cols

    } else {

      stop("Necessary columns for maf file: \n", paste(maf_cols, collapse = " "), "\n")

    }

    if (verbose) cat("-> Keeping only chr 1-23 for maf data...\n")

    maf[, Chromosome := as.character(Chromosome)]
    maf[, Chromosome := gsub(pattern = "chr", replacement = "", Chromosome, ignore.case = TRUE)]
    maf[, Chromosome := gsub(pattern = "X", replacement = "23", Chromosome, ignore.case = TRUE)]
    maf <- maf[Chromosome %in% as.character(1:23), ]

  }


  #-- create temp files to generate input for RunAbsolute

  samples_id <- unique( seg$Sample )

  if (recover) {

    cat("-> recover mode is TRUE, checking samples have been called...\n")

    samples_id <- samples_id[
      !file.exists(
        file.path( cache_dir, paste0(samples_id, ".ABSOLUTE.RData") )
      )
    ]

    if (length(samples_id) != 0) {
      cat("-> ABSOLUTE calling for above samples will be recovered.\n")
    } else {
      cat("-> ABSOLUTE calling has been done. \n")
    }

  }

  seg_filepath <- vector(mode = "character", length = length(samples_id))
  maf_filepath <- vector(mode = "character", length = length(samples_id))

  if (verbose) cat("-> Spliting seg data of samples to different files...\n")

  for (i in seq_along(samples_id)) {

    seg_filepath[i] <- file.path(temp_dir, paste0(samples_id[i], ".seg"))

    data.table::fwrite(
      x = seg[Sample == samples_id[i], ],
      file = seg_filepath[i],
      sep = "\t"
    )

    if ( is.null(maf) ) {

      maf_filepath[i] <- NA_character_

    } else {

      maf_temp <- maf[Tumor_Sample_Barcode == samples_id[i], ]

      if (nrow(maf_temp) == 0) {

        if (verbose) cat("---> This sample has not maf data, skipping...\n")

        maf_filepath[i] <- NA_character_

      } else {

          if (verbose) cat("---> Outputing corresponding maf file...\n")
          data.table::fwrite(
            x = maf_temp,
            file = file.path(temp_dir, paste0(samples_id[i], ".maf")),
            sep = "\t")

      }
    }
  }

  if (verbose) cat("-> Spliting seg data of samples done.\n")

  #-- match options
  platform <- match.arg(platform)
  copy_num_type <- match.arg(copy_num_type)

  if (length(samples_id) != 0) {

    if (verbose) cat("-> Running ABSOLUTE algorithm...", "\n")

    if (is.null(BPPARAM)) BPPARAM <- switch (
      Sys.info()[["sysname"]],
      Linux = BiocParallel::MulticoreParam(progressbar = verbose),
      Windows = BiocParallel::SnowParam(progressbar = verbose)
    )

    BiocParallel::bplapply(
      seq_along(samples_id),
      function(i, ..fn, seg_filepath, maf_filepath,
               samples_id, sigma_p, max_sigma_h, min_ploidy, max_ploidy,
               primary_disease, platform, cache_dir, max_as_seg_count,
               max_non_clonal, max_neg_genome, copy_num_type, min_mut_af,
               error_dir, verbose){

        suppressPackageStartupMessages( loadNamespace("ABSOLUTE") )

        if (is.na(maf_filepath[i])) {
          maf_fn <- NULL
        } else {
          maf_fn <- maf_filepath[i]
        }

        if (verbose) cat("--> Processing sample ", samples_id[i], "...\n")
        ..fn(
          seg_dat_fn = seg_filepath[i], maf_fn = maf_fn,
          sample_name = samples_id[i],
          sigma_p = sigma_p, max_sigma_h = max_sigma_h,
          min_ploidy = min_ploidy, max_ploidy = max_ploidy,
          primary_disease = primary_disease, platform = platform,
          cache_dir = cache_dir, max_as_seg_count = max_as_seg_count,
          max_non_clonal = max_non_clonal, max_neg_genome = max_neg_genome,
          copy_num_type = copy_num_type, min_mut_af = min_mut_af,
          error_dir = error_dir,
          verbose = verbose
        )

      }, ..fn = safe_runAbsolute, seg_filepath = seg_filepath,
      maf_filepath = maf_filepath, samples_id = samples_id,
      sigma_p = sigma_p, max_sigma_h = max_sigma_h,
      min_ploidy = min_ploidy, max_ploidy = max_ploidy,
      primary_disease = primary_disease, platform = platform,
      cache_dir = cache_dir, max_as_seg_count = max_as_seg_count,
      max_non_clonal = max_non_clonal, max_neg_genome = max_neg_genome,
      copy_num_type = copy_num_type, min_mut_af = min_mut_af,
      error_dir = results_dir,
      verbose = verbose, BPPARAM = BPPARAM)

    if (verbose) cat("-> RunAbsolute done. \n\n-> Retrieving results...\n")


    ## This will make problem when iteration computation with same temp dir
    # absolute_files = file.path(cache_dir, grep("RData", dir(cache_dir), value = TRUE))
    absolute_files <- file.path(
      cache_dir, paste0(samples_id, ".ABSOLUTE.RData"))

    cat("--> Files in cache directory:\n")
    print(dir(cache_dir))
    if (verbose) cat("-> Checking result files...\n")

    for (f in absolute_files) {
      if (!file.exists(f)) {
        warning("--> Result file ", f,
                " does not exist, drop it.", immediate. = TRUE)
        absolute_files <- setdiff(absolute_files, f)
      }
    }

    if (length(absolute_files) < 1) {
      stop("No result file to proceed.")
    }

    if (verbose) cat("-> Checked.\n")

    review_dir <- file.path(cache_dir, "review")

    if (dir.exists(review_dir)) {
      if (verbose) cat("-> Removed previous temp review result directory.\n")
      unlink(review_dir, recursive = TRUE)
    }

    if (verbose) cat("-> Running Absolute summarize...\n")
    suppressWarnings(ABSOLUTE::CreateReviewObject(
      obj.name = "DoAbsolute",
      absolute.files = absolute_files,
      indv.results.dir = review_dir,
      copy_num_type = copy_num_type,
      plot.modes = TRUE, verbose = verbose
    ))
    if (verbose) cat("\n-> Absolute summarize done. Prepare auto-reviewing...\n")

    # pp_call_fn = file.path(review_dir, grep("PP-calls_tab.txt", dir(review_dir), value = TRUE))
    # modes_fn = file.path(review_dir, grep("PP-modes.data.RData", dir(review_dir), value = TRUE))
    pp_call_fn <- file.path(review_dir, "DoAbsolute.PP-calls_tab.txt")
    modes_fn <- file.path(review_dir, "DoAbsolute.PP-modes.data.RData")
    suppressWarnings(ABSOLUTE::ExtractReviewedResults(
      reviewed.pp.calls.fn = pp_call_fn,
      analyst.id = "YJ",
      modes.fn = modes_fn,
      out.dir.base = review_dir,
      obj.name = "DoAbsolute",
      copy_num_type = copy_num_type,
      verbose = verbose
    ))
    if (verbose) cat("-> Absolute Auto-reviewing done.\n")

    reviewed_dir <- file.path(review_dir, "reviewed")

    if (verbose) cat("-> Outputing final results...\n")
    if (keep_all_results) {
      cat("--> Choose keeping all results...\n")
    } else {
      cat("--> Choose not keeping all results. Keepping only final results...\n")
    }

    seg_dir <- file.path(results_dir, "seg")
    maf_dir <- file.path(results_dir, "maf")
    dir_create2(seg_dir)
    dir_create2(maf_dir)

    cat("--> Copying DoAbsolute files in review dir to result directory...\n")
    file.copy(
      from = list.files(
        path = reviewed_dir,
        pattern = "DoAbsolute",
        full.names = TRUE),
      to = results_dir
    )

    files_seg <- list.files(
      file.path(reviewed_dir, "SEG_MAF"),
      pattern = "segtab.txt",
      full.names = TRUE
    )
    cat("--> Copying ",
        paste(files_seg, collapse = ", "), "to", seg_dir, "...\n")
    file.copy(from = files_seg, to = seg_dir)


    files_maf <- list.files(
      file.path(reviewed_dir, "SEG_MAF"),
      pattern = "ABS_MAF.txt",
      full.names = TRUE
    )

    cat("--> Copying ",
        paste(files_maf, collapse = ", "), "to", maf_dir, "...\n")
    file.copy(from = files_maf, to = maf_dir)


    if (keep_all_results) {

      summary_res_dir <- file.path(results_dir, "summary_res")
      call_res_dir <- file.path(results_dir, "call_res")
      sample_summaized_dir <- file.path(results_dir, "sample_summaized")

      dir_create2(summary_res_dir)
      dir_create2(call_res_dir)
      dir_create2(sample_summaized_dir)

      files_called_res <- list.files(
        review_dir, pattern = "DoAbsolute",
        all.files = TRUE, full.names = TRUE
      )
      file.copy(from = files_called_res, to = summary_res_dir)


      absolute_plot_files <- file.path(
        cache_dir, paste0(samples_id, ".ABSOLUTE_plot.pdf"))
      absolute_plot_files <- absolute_plot_files[
        file.exists(absolute_plot_files)
      ]
      file.copy(from = absolute_files, to = call_res_dir)
      file.copy(from = absolute_plot_files, to = call_res_dir)



      files_sample_summaized <- list.files(
        file.path(reviewed_dir, "samples"),
        all.files = FALSE, full.names = TRUE
      )
      file.copy(from = files_sample_summaized, to = sample_summaized_dir)

    }

    if (file.exists(file.path(results_dir, "error.log"))) {
      cat("-> Error log info detected, see error.log under result directory for details.\n")
    }

    cat("-> Done.\n")
  }
}


# run_absolute utility functions ----------------------------------------------

safe_runAbsolute <- function(
  seg_dat_fn, maf_fn, sample_name, sigma_p, max_sigma_h,
  min_ploidy, max_ploidy, primary_disease, platform,
  cache_dir, max_as_seg_count, max_non_clonal,
  max_neg_genome, copy_num_type, min_mut_af, error_dir,
  verbose){

  tryCatch(
    {
      suppressWarnings(ABSOLUTE::RunAbsolute(
        seg.dat.fn = seg_dat_fn, maf.fn = maf_fn,
        sample.name = sample_name,
        sigma.p = sigma_p, max.sigma.h = max_sigma_h,
        min.ploidy = min_ploidy, max.ploidy = max_ploidy,
        primary.disease = primary_disease, platform = platform,
        results.dir = cache_dir, max.as.seg.count = max_as_seg_count,
        max.non.clonal = max_non_clonal,
        max.neg.genome = max_neg_genome,
        copy_num_type = copy_num_type,
        min.mut.af = min_mut_af, verbose = verbose
      ))
    },
    error = function(e) {
      cat("----> Error detected, see log for more details.\n")
      sink(file.path(error_dir, "error.log"), append = TRUE)
      cat("Detected error in sample", sample_name, "\n")
      cat("Error message:", e$message, "\n")
      if (grepl("mutations left", e$message)) {
        cat("Try fixing by removing maf file.\n")
        tryCatch(
          {
            suppressWarnings(ABSOLUTE::RunAbsolute(
              seg.dat.fn = seg_dat_fn,  maf.fn = NULL,
              sample.name = sample_name,
              sigma.p = sigma_p, max.sigma.h = max_sigma_h,
              min.ploidy = min_ploidy, max.ploidy = max_ploidy,
              primary.disease = primary_disease, platform = platform,
              results.dir = cache_dir,
              max.as.seg.count = max_as_seg_count,
              max.non.clonal = max_non_clonal,
              max.neg.genome = max_neg_genome,
              copy_num_type = copy_num_type,
              verbose = verbose
            ))
            cat("Fixing successfully!\n")
          },
          error = function(e) {
            cat("Fixing failed. Skipping this sample.\n")
          }
        )
      } else {
        cat("Skipping this sample.\n")
      }
      cat("========\n")
      sink()
    }
  )

}

dir_create2 <- function(x) {
  if (!dir.exists(x)) {
    message("Directory ", x, " does not exists. try to create it")
    if (!dir.create(x, recursive = TRUE)) {
      if (dir.create(x)) {
        stop("Failed creating directory!")
      } else {
        message("Succeeded in creating directory")
      }
    } else {
      message("Succeeded in creating directory")
    }
  }
}
