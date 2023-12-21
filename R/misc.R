#' select columns using tidyselect
#'
#' @param ... \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} One or more
#'   unquoted expressions separated by commas. Variable names can be used as if
#'   they were positions in the data frame, so expressions like `x:y` can be
#'   used to select a range of variables.
#' @param data A data frame or data.table
#' @export
cols_c <- function(..., data) {
  assert_pkg("tidyselect")
  arg_dots <- rlang::enquos(c(...))
  tidyselect::eval_select(arg_dots, data = data)
}
