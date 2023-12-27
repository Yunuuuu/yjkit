`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom rlang .data
NULL

#' format a numeric vector
#'
#' @param x An atomic vector
#' @return A character vector containing the formatted x
#' @export
format_num <- function(x) {
    assert_(x, is.numeric, "numeric")
    data.table::fifelse(
        x < 0.01,
        sprintf("%.2e", x),
        sprintf("%.3f", x)
    )
}

pull <- function(.data, var = -1, name = NULL) {
    var <- tidyselect::vars_pull(names(.data), !!rlang::enquo(var))
    name <- rlang::enquo(name)
    if (rlang::quo_is_null(name)) {
        return(.data[[var]])
    }
    name <- tidyselect::vars_pull(names(.data), !!name)
    rlang::set_names(.data[[var]], nm = .data[[name]])
}
