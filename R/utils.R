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


