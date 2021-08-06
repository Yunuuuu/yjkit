#' Caculate Chromosome-arm-levels copy number variation
#'
#' @param seg_cnv is a GenomicRanges obeject with interger values of CNV column
#'   in \code{\link[S4Vectors:Vector-class]{mcols}(seg_cnv)[[cnv_col]]}
#'   (absolute segment copy number determined by ABSOLUTE algorithm or relative
#'   segment copy number defined with -1 meaning del, 0 meaning neutral and 1
#'   meaning amp) and samples ID in
#'   \code{\link[S4Vectors:Vector-class]{mcols}(seg_cnv)[[sample_id_col]]}
#' @param cnv_col a scalar character gives the column containing CNV values.
#' @param sample_id_col a scalar character gives the column containing sample ID
#' @param ref_cytoband is a GenomicRanges obeject containing the Cytoband
#'   reference, It can be a scalar character \code{"hg19"} or \code{"hg38"}, or
#'   you can provided a self-defined GenomicRanges obeject. \code{"hg38"} is
#'   derived from AnnotationHub by record id \code{"AH53178"} and \code{"hg19"}
#'   is by record id \code{"AH53177"}. Default: \code{"hg38"}.
#' @param filter_centromere Whether to include or exclude segments across
#'   centromere. Default: \code{TRUE}
#' @param cnv_mode is a scala character with values in \code{"rel"} and
#'   \code{"abs"} correspongding to what Shukla, A and Cohen-Sharir have
#'   presented respectively
#' @param threshold the fraction to define Chromosome arm-level aneuploidy
#'   profiling,  Cohen-Sharir uses 0.9 as the cut-off. (\code{cnv_mode}:
#'   "\code{rel}")
#' @param ploidy the ploidy for the comparison to define Chromosome arm-level
#'   aneuploidy profiling. Cohen-Sharir uses background ploidy (\code{cnv_mode}:
#'   "\code{abs}") derived from \code{\link[=run_absolute]{ABSOLUTE}} algorithm.
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @return a tibble containing Chromosome-arm-levels (CNV)
#' @references  \itemize{\item   Cohen-Sharir, Y., McFarland, J.M., Abdusamad,
#'   M. et al. Aneuploidy renders cancer cells vulnerable to mitotic checkpoint
#'   inhibition. Nature 590, 486–491 (2021).
#'   \url{https://doi.org/10.1038/s41586-020-03114-6} \item Shukla, A., Nguyen,
#'   T.H.M., Moka, S.B. et al. Chromosome arm aneuploidies shape tumour
#'   evolution and drug response. Nat Commun 11, 449 (2020).
#'   \url{https://doi.org/10.1038/s41467-020-14286-0}}
#' @examples
#' seg_cnv <- readRDS(system.file("extdata", "run_arm_cnv_example_seg_cnv.rds",
#'                                package = "yjtools"))
#' arm_cnv_res <- run_arm_cnv(seg_cnv, "CNV", "barcode")
#' @export
run_arm_cnv <- function(
  seg_cnv, cnv_col, sample_id_col = NULL,
  ref_cytoband = NULL, filter_centromere = TRUE,
  cnv_mode = c("rel", "abs"), ploidy = NULL,
  threshold = 0.9){

  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)){
    stop("GenomeInfoDb needed for this function to work. Please install it",
         call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)){
    stop("S4Vectors needed for this function to work. Please install it",
         call. = FALSE)
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)){
    stop("GenomicRanges needed for this function to work. Please install it",
         call. = FALSE)
  }

  stopifnot(inherits(seg_cnv, "GenomicRanges"))

  cnv_mode <- match.arg(cnv_mode)

  if (!rlang::is_scalar_character(cnv_col)) stop(
    "cnv_col should be a scalar character", call. = FALSE
  )

  if(!cnv_col %in% colnames(S4Vectors::mcols(seg_cnv))) {
    stop("the column ", cnv_col, " ",
         "is not in S4Vectors::mcols(seg_cnv).",
         call. = FALSE
    )
  }

  if( !is.numeric(S4Vectors::mcols(seg_cnv)[[cnv_col]]) ) {
    stop("the type of CNV column should be a numeric but not a type of ",
         typeof(S4Vectors::mcols(seg_cnv)[[cnv_col]]),
         call. = FALSE)
  }

  stopifnot( is.logical(filter_centromere) )

  if (cnv_mode == "abs") {

    if( !all(S4Vectors::mcols(seg_cnv)[[cnv_col]]  >= 0) ) stop(
      'The CNV values in S4Vectors::mcols(seg_cnv)[[cnv_col]] ',
      'should be above 0 when cnv_mode is "abs"',
      call. = FALSE
    )

    if (is.null(ploidy)) ploidy <- 2L else if( !is.numeric(ploidy) ) {
      stop("The type of ploidy should be a numeric vector but not a ",
           "type of ", typeof(ploidy), ' when cnv_mode is "abs"',
           call. = FALSE
      )
    }

  }

  if (cnv_mode == "rel") {
    if( !all(S4Vectors::mcols(seg_cnv)[[cnv_col]] %in% -1:1) ) stop(
      'The CNV values in S4Vectors::mcols(seg_cnv)[[cnv_col]] ',
      'should be in -1, 0 and 1 when cnv_mode is "rel". ',
      "-1 means Del of copy number, 1 means Amp of copy number ",
      "and 0 means neutral of copy number.",
      call. = FALSE
    )
    if( !(is.numeric(threshold) && length(threshold) == 1) ) stop(
      'threshold should be a scalar numeric when cnv_mode is "abs" ',
      "but not a type of ", typeof(threshold), " ",
      "with length ", length(threshold), ".",
      call. = FALSE
    )
  }


  if ( all(GenomeInfoDb::seqlevelsStyle(seg_cnv) != "UCSC") ){

    message("Mapping seqnames of seg_cnv to UCSC style", appendLF = TRUE)
    GenomeInfoDb::seqlevels(seg_cnv) <- "UCSC"

  }

  if (is.null(ref_cytoband)) ref_cytoband <- run_arm_cnv_ref_cytoband_hg38

  if (rlang::is_scalar_character(ref_cytoband) &&
      ref_cytoband %in% c("hg19", "hg38")){
    ref_cytoband <- rlang::eval_bare(
      rlang::sym(paste0("run_arm_cnv_ref_cytoband_", ref_cytoband))
    )
  } else if (!inherits(ref_cytoband, "GenomicRanges")){
    stop('ref_cytoband should be a scalar character "hg19" or "hg38", ',
         "or a self-defined GenomicRanges obeject.",
         call. = FALSE)
  }
  arm_cytoband <- cytoband_to_arm( ref_cytoband )


  if( !is.null(sample_id_col) ){

    if (!rlang::is_scalar_character(sample_id_col)) stop(
      "sample_id_col should be a scalar character", call. = FALSE
    )

    if(!sample_id_col %in% colnames(S4Vectors::mcols(seg_cnv))) {
      stop("the column ", sample_id_col, " ",
           "is not in S4Vectors::mcols(seg_cnv).",
           call. = FALSE
      )
    }

    sample_id <- unique( S4Vectors::mcols(seg_cnv)[[sample_id_col]] )

    if (cnv_mode == "abs") {

      if( !(length(ploidy) == 1 || length(sample_id) == length(ploidy)) ) stop(
        "the length of ploidy must equal to 1 or the number of unique ",
        "values in the column sample_id_col of seg_cnv.",
        call. = FALSE
      )

      if ( length(ploidy) == 1 ) {
        ploidy <- rep( ploidy, length.out =  length(sample_id) )
        names(ploidy) <- NULL
      }

      if ( is.null(names(ploidy)) ) {
        names(ploidy) <- sample_id
      } else if (!setequal( names(ploidy), sample_id )) {
        stop("If ploidy has names, the names of ploidy should setequal to ",
             "the column sample_id_col of seg_cnv.",
             call. = FALSE)
      }

    }

    if(filter_centromere){

      acen_region <- arm_cytoband[ arm_cytoband$arm == "acen" ]

      # Remove any segment that sligthly overlaps the centromere

      centromere_hits <- GenomicRanges::findOverlaps( seg_cnv, acen_region )

      if ( length(centromere_hits) > 0) seg_cnv <- seg_cnv[
        -S4Vectors::queryHits(centromere_hits)
      ]

    }

    arm_cytoband <- arm_cytoband[
      arm_cytoband$arm %in% c("p", "q")
    ]

    seg_cnv_list <- GenomicRanges::split(
      x = seg_cnv,
      f = S4Vectors::mcols(seg_cnv)[[sample_id_col]],
      drop = TRUE
    )

    arm_cnv_list <- suppressMessages(
      lapply( names(seg_cnv_list), function(i){

        res <- seg_to_arm_cnv(
          seg_cnv = seg_cnv_list[[i]],
          cnv_col = cnv_col,
          arm_cytoband = arm_cytoband,
          cnv_mode = cnv_mode,
          threshold = threshold,
          ploidy = ploidy[[i]]
        )

        res <- dplyr::mutate( res, Sample = !!i, .before = 1)

      })
    )

    return( do.call(rbind, arm_cnv_list) )

  } else {

    if(filter_centromere){

      acen_region <- arm_cytoband[ arm_cytoband$arm == "acen" ]

      # Remove any segment that sligthly overlaps the centromere

      centromere_hits <- GenomicRanges::findOverlaps( seg_cnv, acen_region )

      if ( length(centromere_hits) > 0) seg_cnv <- seg_cnv[
        -S4Vectors::queryHits(centromere_hits)
      ]

    }

    arm_cytoband <- arm_cytoband[
      arm_cytoband$arm %in% c("p", "q")
    ]

    return( seg_to_arm_cnv(
      seg_cnv = seg_cnv, cnv_col = cnv_col,
      arm_cytoband = arm_cytoband,
      cnv_mode = cnv_mode, threshold = threshold,
      ploidy = ploidy[[1]]
    ) )

  }
}

# get_arm_cnv util function -------------------------------------------

cytoband_to_arm <- function(cytoband_genome){

  stopifnot(inherits(cytoband_genome, "GenomicRanges"))

  if ( all( GenomeInfoDb::seqlevelsStyle(cytoband_genome) != "UCSC" ) ){
    message("try to map seqnames of ref_cytoband to UCSC style",
            appendLF = TRUE)

    GenomeInfoDb::seqlevels(cytoband_genome) <- "UCSC"

  }

  cytoband_genome <- GenomeInfoDb::keepSeqlevels(
    cytoband_genome,
    value = stringr::str_subset(
      GenomeInfoDb::seqlevels(cytoband_genome),
      pattern = "^chr((\\d){1,2}|[XY])$"
    ),
    pruning.mode = "tidy"
  )

  S4Vectors::mcols(cytoband_genome)$arm <- factor(
    dplyr::case_when(
      cytoband_genome$gieStain == "acen" ~ "acen",
      TRUE ~ stringr::str_extract( cytoband_genome$name, "[pq]" )
    ),
    levels = c("p", "acen", "q")
  )

  cytoband_genome <- GenomicRanges::split(
    cytoband_genome, GenomeInfoDb::seqnames(cytoband_genome)
  )

  cytoband_genome <- S4Vectors::endoapply(
    cytoband_genome, function(cytoband){

      chro_res <- unlist(GenomicRanges::reduce(
        GenomicRanges::split(cytoband, f = cytoband$arm)
      ), use.names = TRUE)
      S4Vectors::mcols(chro_res)$arm <- factor(
        names(chro_res), levels = c("p", "acen", "q")
      )
      chro_res

    })
  res <- unlist(cytoband_genome, use.names = FALSE)
  names(res) <- NULL

  res
}

seg_to_arm_cnv <- function(
  seg_cnv, cnv_col, arm_cytoband, cnv_mode = "rel", threshold, ploidy
){

  stopifnot( inherits(seg_cnv, "GenomicRanges") )

  stopifnot( inherits(arm_cytoband, "GenomicRanges") )

  S4Vectors::mcols(arm_cytoband)$arm_width <-
    GenomicRanges::width( arm_cytoband )

  overlap_hits <- GenomicRanges::findOverlaps(
    seg_cnv, arm_cytoband, type = "any"
  )

  intersect_region <- GenomicRanges::pintersect(
    seg_cnv[S4Vectors::queryHits(overlap_hits)],
    arm_cytoband[S4Vectors::subjectHits(overlap_hits)],
    drop.nohit.ranges = FALSE,
    ignore.strand = FALSE,
    strict.strand = FALSE
  )

  res <- tibble::tibble(
    seqnames = as.factor( GenomicRanges::seqnames(intersect_region) ),
    width = GenomicRanges::width(intersect_region),
    CNV = abs( S4Vectors::mcols( intersect_region )[[cnv_col]] )
  ) %>%
    dplyr::mutate(

      arm = !!S4Vectors::mcols(
        arm_cytoband[ S4Vectors::subjectHits(overlap_hits) ]
      )[["arm"]],

      arm_width = !!S4Vectors::mcols(
        arm_cytoband[ S4Vectors::subjectHits(overlap_hits) ]
      )[["arm_width"]],

      seg_length_frac = .data$width / .data$arm_width

    ) %>%
    dplyr::group_by(.data$seqnames, .data$arm)

  if (cnv_mode == "rel") res <- res %>%
    dplyr::summarize(
      arm_cnv = as.integer(
        sum(.data$CNV * .data$seg_length_frac, na.rm = TRUE) > !!threshold
      )
    )

  if (cnv_mode == "abs") res <- res %>%
    dplyr::summarize(
      arm_cnv = as.integer(sign(
        round(matrixStats::weightedMedian(
          .data$CNV, w = .data$width
        ), digits = 0) - round(!!ploidy, digits = 0)
      ))
    )

  dplyr::ungroup(res)

}
