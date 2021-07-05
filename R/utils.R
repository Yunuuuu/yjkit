#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' format a numeric vector
#'
#' @param x a numeric vector
#' @return a character vector containing the formatted x
format_num <- function(x){

  stopifnot(is.numeric(x))
  dplyr::if_else(
    x < 0.01,
    sprintf("%1.2e", x),
    sprintf("%.3f", x)
  )

}


#' select columns using tidyselect
#'
#' @param ... <\code{\link[dplyr:dplyr_tidy_select]{tidy-select}}> One or more
#'   unquoted expressions separated by commas. Variable names can be used as if
#'   they were positions in the data frame, so expressions like `x:y` can be
#'   used to select a range of variables.
#' @param data 	A data frame or data.table
#' @import tidyselect
#' @export
cols_c <- function(..., data){

  arg_dots <- rlang::expr(c(...))
  tidyselect::eval_select(arg_dots, data = data)

}
