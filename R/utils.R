#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom rlang .data
NULL

#' format a numeric vector
#'
#' @param x a atomic vector
#' @return a character vector containing the formatted x
format_num <- function(x){

  if(is.numeric(x)){
    x <- dplyr::if_else(
      x < 0.01,
      sprintf("%1.2e", x),
      sprintf("%.3f", x)
    )
  }

  if (is.factor(x)) x <- as.character(x)

  x
}


