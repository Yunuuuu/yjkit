#' select columns using tidyselect
#'
#' @param ... \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} One or more
#'   unquoted expressions separated by commas. Variable names can be used as if
#'   they were positions in the data frame, so expressions like `x:y` can be
#'   used to select a range of variables.
#' @param data A data frame or data.table
#' @import tidyselect
#' @export
cols_c <- function(..., data){

  arg_dots <- rlang::expr(c(...))
  tidyselect::eval_select(arg_dots, data = data)

}


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

  if (check_logged(expr_data)) {

    print("log2 transformation wasn't needed")
    return(data)

  } else {

    print("Doing log2 transformation")
    expr_data <- log2(expr_data + 1)

  }

  if ( inherits(data, "ExpressionSet") ) {

    Biobase::exprs(data) <- expr_data
    return(data)

  }

  if ( inherits(data, "matrix") ) return(expr_data)

}



# check whether vector is log transformation ------------------------------
# a scalar logical value, \code{TRUE} means logged, and \code{FALSE} means not.
check_logged <- function(x){

  qx <- as.numeric(
    stats::quantile(x, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  )
  not_log <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

  !not_log

}
