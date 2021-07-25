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


#' Install or update CRAN, Bioconductor and GitHub packages
#'
#' Install or update CRAN, Bioconductor and GitHub packages from source and
#' compiling packages
#' @inheritParams BiocManager::install
#' @inherit BiocManager::install details
#' @inherit BiocManager::install return
#' @inherit BiocManager::install seealso
#' @export
install_pkgs <- function(pkgs = character(), ..., type = "source",
                    INSTALL_opts = "--byte-compile",
                    site_repository  = character(),
                    update = TRUE, ask = TRUE,
                    checkBuilt = FALSE){

  if (!requireNamespace("BiocManager", quietly = TRUE)){
    stop("BiocManager needed for this function to work. Please install it",
         call. = FALSE)
  }

  BiocManager::install(
    pkgs = pkgs, type = type,
    INSTALL_opts = INSTALL_opts, ...,
    site_repository = site_repository,
    update = update, ask = ask,
    checkBuilt = checkBuilt
  )

  pkgs

}
