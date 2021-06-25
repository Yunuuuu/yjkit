#' auto log2 transform for GEO datasets
#'
#' @param data gene expression data, can be a \code{ExpressionSet} object  or
#'   \code{matrix} object.
#' @details Automatically check whether \code{data} has experientced logarithm
#'   transformation, if not, applying a log2 transformation. The test
#'   methodology for logarithm transformation is based on
#'   \href{https://www.ncbi.nlm.nih.gov/geo/geo2r/}{GEO2R}
#' @references \href{https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE1122}{GEO
#'   analysis}
#' @return a \code{matrix} or \code{ExpressionSet} object of logarithm
#'   transformed if it is unlogged
#'   \code{class(data_set)}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @export
auto_log <- function(data) {

  if (!requireNamespace("Biobase", quietly = TRUE)){
    stop("Biobase needed for this function to work. Please install it",
         call. = FALSE)
  }

  if (!(inherits(data, "ExpressionSet") || inherits(data, "matrix"))) {
    stop("data must be a class of matrix or ExpressionSet")
  }

  if ( inherits(data, "ExpressionSet") ) expr_data <- Biobase::exprs(data)
  if ( inherits(data, "matrix") ) expr_data <- data

  qx <- as.numeric(
    stats::quantile(expr_data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  )
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

  if (LogC) {

    expr_data <- log2(expr_data + 1)

    print("log2 transformation is done")

  } else {

    print("log2 transformation wasn't needed")
    return(data)

  }

  if ( inherits(data, "ExpressionSet") ) {

    Biobase::exprs(data) <- expr_data
    return(data)

  }

  if ( inherits(data, "matrix") ) return(expr_data)

}


