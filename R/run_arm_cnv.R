#' Caculate Chromosome-arm-levels copy number variation
#'
#' @param seg_cnv is a GenomicRanges obeject with interger values of CNV column
#'   in \code{\link[S4Vectors:Vector-class]{mcols}(seg_cnv)[[cnv_col]]}
#'   (absolute segment copy number determined by ABSOLUTE algorithm or relative
#'   segment copy number defined with -1 meaning del, 0 meaning neutral and 1
#'   meaning amp) and samples ID in
#'   \code{\link[S4Vectors:Vector-class]{mcols}(seg_cnv)[[sample_id_col]]}
#' @param cnv_col a scalar character gives the column containing CNV
#' @param sample_id_col a scalar character gives the column containing sample ID
#' @param ref_cytoband is a GenomicRanges obeject containing the Cytoband
#'   reference, It can be a scalar character "hg19" or "hg38", or you can
#'   provided a self-defined GenomicRanges obeject. Default: "hg38"
#' @param filter_centromere Whether to include or exclude segments across
#' centromere. Default: \code{TRUE}
#' @param cnv_mode is a scala character with values in "rel" and "abs"
#'   correspongding to what Shukla, A and Cohen-Sharir have presented
#'   respectively
#' @param threshold the fraction to define Chromosome arm-level aneuploidy
#'   profiling,  Cohen-Sharir uses 0.9 as the cut-off. (\code{cnv_mode
#'   \strong{rel}})
#' @param ploidy the ploidy for the comparison to define Chromosome arm-level
#'   aneuploidy profiling. Cohen-Sharir uses background ploidy (\code{cnv_mode
#'   \strong{abs}})
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
#' package = "yjkit"))
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

  if(
    is.character(cnv_col) &&
    !cnv_col %in% colnames( S4Vectors::mcols(seg_cnv) )
  ) {
    stop(
      paste("the character", cnv_col,
            "is not in S4Vectors::mcols(seg_cnv)")
    )
  }

  if( !is.numeric(S4Vectors::mcols(seg_cnv)[[cnv_col]]) ) {
    stop("the type of CNV column must be numeric")
  }

  stopifnot( is.logical(filter_centromere) )

  if (cnv_mode == "abs") {
    stopifnot( all(S4Vectors::mcols(seg_cnv)[[cnv_col]]  >= 0) )
    stopifnot( is.null(ploidy) || is.numeric(ploidy) )
  }
  if (cnv_mode == "rel") {
    stopifnot( all(S4Vectors::mcols(seg_cnv)[[cnv_col]] %in% -1:1) )
    stopifnot( is.numeric(threshold) && length(threshold) == 1 )
  }


  if ( all(GenomeInfoDb::seqlevelsStyle(seg_cnv) != "UCSC") ){
    message("try to map seqnames of seg_cnv to UCSC style")
    seg_cnv <- GenomeInfoDb::renameSeqlevels(
      seg_cnv, stats::na.omit(
        GenomeInfoDb::mapSeqlevels(
          GenomeInfoDb::seqlevels(seg_cnv), "UCSC") )
    )
  }

  if (is.null(ploidy)) ploidy <- 2L else {
    stopifnot( is.numeric(ploidy) )
  }

  if (is.null(ref_cytoband)) ref_cytoband <- run_arm_cnv_ref_cytoband_hg38

  if (is.character(ref_cytoband) && (length(ref_cytoband) == 1) && ref_cytoband %in% c("hg19", "hg38")){
    ref_cytoband <- rlang::eval_bare(
      rlang::sym(paste0("run_arm_cnv_ref_cytoband_", ref_cytoband))
    )
  }
  arm_cytoband <- cytoband_to_arm( ref_cytoband )


  if( !is.null(sample_id_col) ){

    stopifnot(
      ( is.character(sample_id_col) || is.numeric(sample_id_col) ) &&
        length(sample_id_col) == 1
    )

    if(
      is.character(sample_id_col) &&
      !sample_id_col %in% colnames( S4Vectors::mcols(seg_cnv) )
    ) {
      stop(
        paste("the character", sample_id_col,
              "is not in S4Vectors::mcols(seg_cnv)")
      )
    }

    sample_len <- length(unique( S4Vectors::mcols(seg_cnv)[[sample_id_col]] ))

    if( !(length(ploidy) == 1 || sample_len == length(ploidy)) ) stop(
      paste("the length of ploidy must equal to 1 or the number of samples")
    )

    if ( length(ploidy) == 1 ) {
      ploidy <- rep( ploidy, length.out = sample_len )
      names(ploidy) <- NULL
    }

    if ( is.null(names(ploidy)) ) names(ploidy) <- unique(
      S4Vectors::mcols(seg_cnv)[[sample_id_col]]
    )

    if (!setequal(
      names(ploidy),
      unique(S4Vectors::mcols(seg_cnv)[[sample_id_col]])
    )) {
      stop("the names of ploidy must equal to all the ID in sample_id_col")
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
    message("try to map seqnames to UCSC style")
    cytoband_genome <- GenomeInfoDb::renameSeqlevels(
      cytoband_genome, stats::na.omit(
        GenomeInfoDb::mapSeqlevels(
          GenomeInfoDb::seqlevels(cytoband_genome), "UCSC") )
    )
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

      seg_length_frac = width / arm_width

    ) %>%
    dplyr::group_by(seqnames, arm) %>%
    dplyr::summarize(
      arm_cnv = dplyr::case_when(
        !!cnv_mode == "rel" ~ as.integer(
          sum(CNV * seg_length_frac, na.rm = TRUE) > !!threshold
        ),
        !!cnv_mode == "abs" ~ as.integer(sign(
          round(matrixStats::weightedMedian(
            CNV, w = width
          ), digits = 0) - round(!!ploidy, digits = 0)
        ))
      )
    ) %>%
    dplyr::ungroup()
  res
}
