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
#' @param repos (Optional) character(1) vector representing an additional
#'   repository in which to look for packages to install. This repository will
#'   be prepended to the default repositories (which you can see with
#'   \code{BiocManager::\link[BiocManager]{repositories}}.
#' @param type character, indicating the type of package to download and
#'   install. Will be "\code{source}" except on Windows and some macOS builds:
#'   see the section on ‘Binary packages’ for those.
#' @param INSTALL_opts an optional character vector of additional option(s) to
#'   be passed to \code{R CMD INSTALL} for a source package install. E.g.,
#'   \code{c("--html", "--no-multiarch", "--no-test-load")}. \cr Can also be a
#'   named list of character vectors to be used as additional options, with
#'   names the respective package names.
#' @inherit BiocManager::install details
#' @return  return the pkgs argument invisibly.
#' @seealso \code{\link[BiocManager]{install}} and
#'   \code{\link[utils]{install.packages}}
#' @export
install_pkgs <- function(pkgs = character(), ..., type = "source",
                    INSTALL_opts = "--byte-compile",
                    repos  = character(),
                    update = TRUE, ask = TRUE,
                    checkBuilt = FALSE){

  if (!requireNamespace("BiocManager", quietly = TRUE)){
    stop("BiocManager needed for this function to work. Please install it",
         call. = FALSE)
  }

  BiocManager::install(
    pkgs = pkgs, type = type,
    INSTALL_opts = INSTALL_opts, ...,
    site_repository = repos,
    update = update, ask = ask,
    checkBuilt = checkBuilt
  )

  invisible(pkgs)

}
