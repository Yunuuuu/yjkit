
#' Single-sample GSEA
#'
#' Project each sample within a data set onto a space of gene set enrichment
#' scores using the ssGSEA projection methodology. According to GenePattern
#' (modified) and GSVA
#'
#' @param data_set gene expression data, can be \code{ExpressionSet} object,
#'   \code{SummarizedExperiment}, or \code{matrix} object. when data_set is a
#'   \code{SummarizedExperiment} object, the \code{assay} argument will be used
#'   to extracted the assay with exression matrix (see
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}), if missing
#'   \code{assay}, the first \code{assay} will be used. For ssGSEA, normalizd
#'   exression values with gene length adjusted were needed.
#' @param gene_set_list gene sets can be provided as \code{GeneSetCollection}
#'   object or \code{list} object
#' @param assay used to extracted the assay with exression matrix (see
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}), if missing
#'   \code{assay}, the first \code{assay} will be used.
#' @param NES whether to calculate normalized enrichment scores
#' @param num_perm the number of permutations to calculate NES
#' @param gene_set_selection  list of gene set names on which to project the
#'   input expression data.  Alternatively, this field may be set to ALL,
#'   indicating that the input expression dataset is to be projected to all gene
#'   sets defined in the specified gene set database(s). (default: \code{ALL})
#' @param min_sz minimal size of each geneSet for analyzing
#' @param max_sz maximal size of genes annotated for testing
#' @param sample_norm_type Normalization method applied to expression data.
#'   Supported methods are rank, log.rank, and log.  (Default: rank)
#' @param weight Exponential weight employed in calculation of enrichment
#'   scores.  The default value of \code{0.75} was selected after extensive
#'   testing.  The module authors strongly recommend against changing from
#'   default. (Default: \code{0.75})
#' @param min_overlap min overlap required between genes in gene set and genes
#'   in input (feature dataset) file in order to include that gene set in data
#'   set projection
#' @param BPPARAM see \code{\link[BiocParallel]{bpparam}}. Default: \code{NULL}
#'   means \code{\cr switch (Sys.info()[["sysname"]], \cr Linux =
#'   BiocParallel::MulticoreParam(), \cr Windows = BiocParallel::SnowParam()) }
#' @param verbose	if TRUE, print extra infos. Default: \code{TRUE}
#' @return a \code{list} or \code{ExpressionSet} or \code{SummarizedExperiment}
#'   object of projection results for each gene set and each sample based on
#'   \code{class(data_set)}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @references \itemize{\item Subramanian A, Tamayo P, Mootha VK, Mukherjee S,
#'   Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES,
#'   Mesirov JP. Gene set enrichment analysis: A knowledge-based approach for
#'   interpreting genome-wide expression profiles. PNAS.
#'   2005;102(43):15545-15550. \url{https://doi.org/10.1073/pnas.0506580102}
#'   \item Barbie, D., Tamayo, P., Boehm, J. et al. Systematic RNA interference
#'   reveals that oncogenic KRAS-driven cancers require TBK1. Nature 462,
#'   108–112 (2009). \url{https://doi.org/10.1038/nature08460}}
#' @export
run_ssgsea <- function(

  # gene expression data
  data_set,

  # list of gene set
  gene_set_list,

  assay,

  NES = TRUE, num_perm = 200,

  # "ALL" or list with names of gene sets
  gene_set_selection  = "ALL",
  min_sz = 1, max_sz = Inf,

  # normalization method applied to input feature data:
  # "none", "rank", "log" or "log.rank"
  sample_norm_type = c("rank", "log", "log.rank", "none"),

  # exponential weight applied to ranking in calculation of
  # enrichment score
  weight = 0.75,

  # min overlap required between genes in gene set and genes in input
  # (feature dataset) file in order to include that gene set in data set
  # projection
  min_overlap = 10,

  BPPARAM = NULL,

  verbose = TRUE) {

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)){
    stop("SummarizedExperiment needed for this function to work. Please install it",
         call. = FALSE)
  }

  if (!requireNamespace("GSEABase", quietly = TRUE)){
    stop("GSEABase needed for this function to work. Please install it",
         call. = FALSE)
  }
  if (!requireNamespace("Biobase", quietly = TRUE)){
    stop("Biobase needed for this function to work. Please install it",
         call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)){
    stop("S4Vectors needed for this function to work. Please install it",
         call. = FALSE)
  }
  if (!requireNamespace("BiocParallel", quietly = TRUE)){
    stop("BiocParallel needed for this function to work. Please install it",
         call. = FALSE)
  }

  if (!requireNamespace("GSVA", quietly = TRUE)){
    stop("GSVA needed for this function to work. Please install it",
         call. = FALSE)
  }

  ## arguments-------------------------------------------

  sample_norm_type <- match.arg(sample_norm_type)

  ## validate input parameters--------------------------

  if (sample_norm_type != "none" && sample_norm_type != "rank" && sample_norm_type != "log" && sample_norm_type != "log.rank") {
    stop("invalid value for sample_norm_type argument: ", sample_norm_type)
  }

  stopifnot(inherits(data_set, c("ExpressionSet", "SummarizedExperiment", "matrix")))
  stopifnot(inherits(gene_set_list, c("GeneSetCollection", "list")))


  ## extract gene expression matrix-----------------------------------

  if (inherits(data_set, "ExpressionSet")) expr <- Biobase::exprs(data_set)

  if (inherits(data_set, "SummarizedExperiment")) {

    if (length(SummarizedExperiment::assays(data_set)) == 0L)
      stop("The input SummarizedExperiment object has no assay data.")

    if (missing(assay)) {

      assay <- names(SummarizedExperiment::assays(data_set))[[1]]

    } else {

      if (!is.character(assay))
        stop("The 'assay' argument must contain a character string.")

      assay <- assay[[1]]

      if (!assay %in% names(SummarizedExperiment::assays(data_set)))
        stop(sprintf("Assay %s not found in the input SummarizedExperiment object.", assay))

    }

    expr <- SummarizedExperiment::assays(data_set)[[assay]]

  }

  if ( inherits(data_set, "matrix") ) expr <- data_set

  if (nrow(expr) < 2)
    stop("Less than two genes in the input ExpressionSet object\n")


  ## keep gene and sample names -------------------------------------------

  keep_gene_name <- rownames(expr)
  keep_sample_name <- colnames(expr)


  if ( (inherits(data_set, "ExpressionSet") || inherits(data_set, "SummarizedExperiment")) && inherits(gene_set_list, "GeneSetCollection") ){

    ## map gene identifiers of the gene sets to the features in the chip
    if (inherits(data_set, "ExpressionSet"))
      annotpkg <- Biobase::annotation(data_set)
    if (inherits(data_set, "SummarizedExperiment"))
      annotpkg <- S4Vectors::metadata(data_set)$annotation

    if (!is.null(annotpkg) && length(annotpkg) > 0 && is.character(annotpkg) && annotpkg != "") {

      if (!annotpkg %in% utils::installed.packages())
        stop(sprintf("Please install the nnotation package %s", annotpkg))

      if (verbose) cat("Mapping identifiers between gene sets and feature names\n")

      ## Biobase::annotation() is necessary to disambiguate from the
      ## 'annotation' argument
      mapped_gs_list<- GSEABase::mapIdentifiers(
        gene_set_list,
        GSEABase::AnnoOrEntrezIdentifier(annotpkg)
      )
      mapped_gs_list <- GSEABase::geneIds(mapped_gs_list)

    } else {

      mapped_gs_list <- GSEABase::geneIds(gene_set_list)
      if (verbose) {
        cat(
          "No annotation package name available in the input 'ExpressionSet' object 'expr'.",
          "Attempting to directly match identifiers in 'expr' to gene sets.",
          sep="\n")
      }

    }
  } else if ( inherits(gene_set_list, "GeneSetCollection") ) {
    mapped_gs_list <- GSEABase::geneIds(gene_set_list)
  } else {
    mapped_gs_list <- gene_set_list
  }

  ## remove gene sets from the analysis for which no features are available
  ## and meet the minimum and maximum gene-set size specified by the user

  gene_set_list <- GSVA::filterGeneSets(
    mapped_gs_list,
    min.sz = max(1, min_sz, na.rm = TRUE),
    max.sz = max_sz
  )


  #############################################
  # Sample normalization
  #############################################

  if (sample_norm_type == "none") {

    print(paste("No normalization to be made"))

  } else if (sample_norm_type == "rank") {

    expr <- apply(expr, 2, function(x) {
      as.integer(rank(x, ties.method = "average"))
    })
    expr <- 10000*expr/nrow(expr)

  } else if (sample_norm_type == "log.rank") {

    expr <- apply(expr, 2, function(x) {
      as.integer(rank(x, ties.method = "average"))
    })
    expr <- log(10000*expr/nrow(expr) + exp(1))

  } else if (sample_norm_type == "log") {

    expr[expr < 1] <- 1
    expr <- log(expr + exp(1))

  }

  #############################################
  # Select desired gene sets
  #############################################

  N.gs <- length(gene_set_list)
  if (gene_set_selection[1] != "ALL") {
    locs <- match(gene_set_selection, names(gene_set_list))
    N.gs <- sum(!is.na(locs))
    if (N.gs == 0) {
      stop("No matches with gene_set_selection")
    }
    gene_set_list <- gene_set_list[locs[!is.na(locs)]]
  }


  # eliminate gene sets for which there was insufficient overlap

  print(paste("N.gs before overlap prunning:", N.gs))

  keep_gene_set <- lapply(gene_set_list, function(x, y){
    intersect(x, y)
  }, y = keep_gene_name)
  keep_gene_set <- lapply(
    keep_gene_set, function(gene_set, min){
      length(gene_set) > min
    }, min = min_overlap
  )
  keep_gene_set <- unlist(keep_gene_set)

  N.gs <- sum(keep_gene_set)

  print(paste("N.gs after overlap prunning:", N.gs))

  gene_set_list <- gene_set_list[keep_gene_set]


  ## recover gene and sample names -------------------------------------

  rownames(expr) <- keep_gene_name
  colnames(expr) <- keep_sample_name

  ## run ssGSEA analysis --------------------------------------------------
  # projecting expression data to gene_set_list
  # data.matrix containing gene expression data
  # exponential weight applied to ranking in calculation of enrichment score
  #################################################

  if (is.null(BPPARAM)) BPPARAM <- switch (
    Sys.info()[["sysname"]],
    Linux = BiocParallel::MulticoreParam(progressbar = verbose),
    Windows = BiocParallel::SnowParam(progressbar = verbose)
  )

  ssgsea_es <- project_to_geneset(
    data_matrix = expr,
    gene_set_list = gene_set_list,
    weight = weight,
    NES = NES, num_perm = num_perm,
    BPPARAM = BPPARAM,
    verbose = verbose
  )

  if ( inherits(data_set, "ExpressionSet") ) {

    rval <- ssgsea_es
    names(rval)[names(rval) == "ES"] <- "exprs"
    rval <- methods::new(
      "ExpressionSet",
      assayData = rval,
      phenoData = Biobase::phenoData(data_set),
      experimentData = Biobase::experimentData(data_set),
      annotation = ""
    )
  }

  if ( inherits(data_set, "SummarizedExperiment") ) {

    rval <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(ssgsea_es),
      colData = SummarizedExperiment::colData(data_set),
      metadata = S4Vectors::metadata(data_set)
    )

  }

  if ( inherits(data_set, "matrix" ) ) {

    rval <- ssgsea_es

  }
  return(rval)
} # ssGSEA according to GenePattern and GSVA


# utility functions -------------------------------------------------------

## optimized version of the function .rndWalk by Alexey Sergushichev
## https://github.com/rcastelo/GSVA/pull/15
## based on his paper https://doi.org/10.1101/060012

gsea <- function(gene_list, gene_set, weight) {

  n <- length(gene_list)
  k <- length(gene_set)
  idx <- sort.int(fastmatch::fmatch(gene_set, names(gene_list)))
  gene_list <- abs(gene_list)^weight
  stepCDFinGeneSet2 <-
    sum(gene_list[idx] * (n - idx + 1)) / sum(gene_list[idx])

  stepCDFoutGeneSet2 <- (n * (n + 1) / 2 - sum(n - idx + 1)) / (n - k)

  walkStat <- stepCDFinGeneSet2 - stepCDFoutGeneSet2

  walkStat

}

perm_gene_set <- function(gene_list) {
  perm.idx <- sample.int(length(gene_list), replace = FALSE)
  perm_gene_set <- gene_list
  names(perm_gene_set) <- names(gene_list)[perm.idx]
  return(perm_gene_set)
}

perm_gsea_es <- function(gene_list, gene_set, weight) {

  gene_list <- perm_gene_set(gene_list)
  res <- gsea(
    gene_list = gene_list,
    gene_set = gene_set,
    weight = weight
  )
  return(res)

}

project_to_geneset <- function(data_matrix, gene_set_list, weight, NES, num_perm, BPPARAM, verbose){

  gene_name <- rownames(data_matrix)
  message("Running ssGSEA.............")
  es_res <- BiocParallel::bplapply(
    seq_len(ncol(data_matrix)),
    function(j, data_matrix, gene_set_list, weight, NES,
             num_perm, verbose) {

      gene_list <- data_matrix[, j, drop = TRUE]
      gene_list <- (gene_list-mean(gene_list))/stats::sd(gene_list)
      names(gene_list) <- gene_name
      gene_list <- sort(gene_list, decreasing = TRUE)

      if (verbose) {
        message_nes <- ""
        if (NES) message_nes <- ", NES and p-values"
        message(
          "calculating ES", message_nes, " for ",
          j, "th sample ", colnames(data_matrix)[j]
        )
      }

      enrich_score <- lapply(gene_set_list, function(gene_set){

        es <- gsea(
          gene_list = gene_list,
          gene_set = gene_set,
          weight = weight
        )

        if (NES) {
          perm_scores <- lapply(1:num_perm, function(i) {
            perm_gsea_es(
              gene_list = gene_list,
              gene_set = gene_set,
              weight = weight
            )
          })

          perm_scores <- unlist(perm_scores)
          pos_m <- mean(perm_scores[perm_scores > 0])
          neg_m <- abs(mean(perm_scores[perm_scores < 0]))


          if (sign(es) == 1 || sign(es) == 0) m <- pos_m
          if (sign(es) == -1) m <- neg_m

          NES <- es/m

          if( is.na(NES) ) {
            pvals <- NA
          } else if ( NES >= 0 ) {
            pvals <- (sum(perm_scores >= NES) + 1) / (sum(perm_scores >= 0)+1)
          } else { # NES < 0
            pvals <- (sum(perm_scores <= NES) + 1) / (sum(perm_scores < 0)+1)
          }

          return(list(ES = es, NES = NES, P.values = pvals))
        }

        return(list(ES = es))

      })

      return(enrich_score)

    },
    data_matrix = data_matrix,
    gene_set_list = gene_set_list,
    weight = weight, NES = NES, num_perm = num_perm,
    BPPARAM = BPPARAM)

  names(es_res) <- colnames(data_matrix)

  ES <- lapply(es_res, function(smp) {

    sample_es_vector <- lapply(smp, function(gene_set) gene_set$ES)
    return(unlist(sample_es_vector, use.names = TRUE))

  })
  ES <- do.call("cbind", ES)
  res <- list(ES = ES)

  if (NES) {
    NES <- lapply(es_res, function(smp) {
      sample_nes_vector <- lapply(smp, function(gene_set) gene_set$NES)
      return(unlist(sample_nes_vector, use.names = TRUE))
    })
    NES <- do.call("cbind", NES)

    P.values <- lapply(es_res, function(smp) {
      sample_pva_vector <- lapply(smp, function(gene_set) {
        return(gene_set$P.values)
      })
      return(unlist(sample_pva_vector, use.names = TRUE))
    })
    P.values <- do.call("cbind", P.values)

    res$NES <- NES
    res$P.values <- P.values
  }
  return(res)

}
