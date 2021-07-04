#' Between group test
#'
#' Function to calculate statistics and P-value for both parametric and
#' nonparametric tests.
#'
#' @param data a data.frame object containing the variables in the x and y
#' @param x,y quoted arguments, can be a list, if x or y is a list object, each
#'   of x and each of y will be paired to implement test
#' @param type if \code{type = "nonparametric"},
#'   \code{\link[stats]{wilcox.test}} or \code{\link[stats]{kruskal.test}} is
#'   used. \cr if \code{type = "parametric"}, \code{\link[stats]{t.test}} or
#'   \code{\link[stats]{aov}} is used.\cr which method is used depends on the
#'   unique values of \code{x} in \code{data}
#' @param ... other arguments passed to statistical test function,
#'   \code{\link[stats]{wilcox.test}}, \code{\link[stats]{kruskal.test}},
#'   \code{\link[stats]{t.test}} and \code{\link[stats]{aov}}
#' @section Quasiquotation:
#'   \code{x} and \code{y} in \code{stat_between_test()}
#'   are \code{\link[rlang:nse-defuse]{quoted arguments}} This means that its
#'   inputs are quoted to be evaluated in the context of the data. This makes it
#'   easy to work with variables from the data frame because you can name those
#'   directly. The flip side is that you have to use
#'   \code{\link[rlang:nse-force]{quasiquotation}} to program with
#'   \code{stat_between_test()}. See a tidy evaluation tutorial such as the
#'   \href{https://dplyr.tidyverse.org/articles/programming.html}{dplyr
#'   programming vignette} to learn more about these techniques.
#'
#' @return a tibble of the test results where each row match each pair of x and
#'   y. (see \code{\link[generics:tidy]{tidy}})
#' @details if \code{type = "nonparametric"}, \code{\link[stats]{wilcox.test}}
#'   or \code{\link[stats]{kruskal.test}} is used. \cr if \code{type =
#'   "parametric"}, \code{\link[stats]{t.test}} or \code{\link[stats]{aov}} is
#'   used.\cr which method is used depends on the unique values of \code{x term}
#'   in \code{data}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#'   stat_between_test(mtcars, vs, mpg, type = "p")
#'   stat_between_test(mtcars, vs, mpg, type = "n")
#'   stat_between_test(mtcars, cyl, mpg, type = "p")
#'   stat_between_test(mtcars, cyl, mpg, type = "n")
#' @export
stat_between_test <- function(data, x, y,
                              type = c("nonparametric", "parametric"),
                              ...){

  stopifnot(inherits(data, "data.frame"))
  type <- match.arg(type)

  quo_list_x <- rlang::enquos(x)
  label_list_x <- rlang::enexprs(x)
  quo_list_y <- rlang::enquos(y)
  label_list_y <- rlang::enexprs(y)

  purrr::map_dfr(seq_along(quo_list_y), function(y_i){

    purrr::map_dfr(seq_along(quo_list_x), function(x_i){

      stat_between_test_helper(data = data,
                             x = !!quo_list_x[[x_i]],
                             y = !!quo_list_y[[y_i]],
                             x_label = deparse1(label_list_x[[x_i]]),
                             y_label = deparse1(label_list_y[[y_i]]),
                             type = type,)

    })
  })


}

#' Correlation analysis
#'
#' Function to calculate Correlation coefficient and P-value between the columns
#' of x and the columns of y, using one of Pearson's product moment correlation
#' coefficient, Kendall's tau or Spearman's rho.
#'
#' @param x,y data.frame object. x and y must have the same \code{nrow}
#' @param use an optional character string giving a method for computing
#'   covariances in the presence of missing values. This must be (an
#'   abbreviation of) one of the strings \code{"everything"}, \code{"all.obs"},
#'   \code{"complete.obs"}, \code{"na.or.complete"}, or
#'   \code{"pairwise.complete.obs"}. Default: \code{"pairwise.complete.obs"}.
#'   See \code{\link[stats]{cor}}
#' @param method a character string indicating which correlation coefficient (or
#'   covariance) is to be computed. One of \code{"pearson"}, \code{"kendall"},
#'   or \code{"spearman"}, can be abbreviated. Default: \code{"spearman"}
#' @param alternative indicates the alternative hypothesis and must be one of
#'   \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just
#'   the initial letter. \code{"greater"} corresponds to positive association,
#'   \code{"less"} to negative association.
#' @param exact a logical indicating whether an exact p-value should be
#'   computed. Used for Kendall's tau and Spearman's rho. See
#'   \code{\link[stats]{cor.test}}. Default: \code{TRUE}.
#' @param conf_level confidence level for the returned confidence interval.
#'   Currently only used for the Pearson product moment correlation coefficient
#'   if there are at least 4 complete pairs of observations.
#' @param ... Other parameters for \code{\link[stats]{cor.test}}
#' @return a tibble of the Correlation results. (see
#'   \code{\link[stats]{cor.test}})
#' @details See \code{\link[stats]{cor.test}}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#'   stat_cor_test(mtcars)
#' @export
stat_cor_test <- function(x, y = NULL,
                          use = c("pairwise.complete.obs",
                                  "everything", "all.obs", "complete.obs",
                                  "na.or.complete"),
                          method = c("spearman", "pearson", "kendall"),
                          alternative = c("two.sided", "less", "greater"),
                          exact = TRUE,
                          conf_level = 0.95, ...){

  if(is.null(y)) y <- x

  stopifnot(inherits(x, "data.frame"))
  stopifnot(inherits(y, "data.frame"))

  use <- match.arg(use)
  method <- match.arg(method)

  cor_res <- purrr::map(x, function(x_sg){

    purrr::map(y, function(y_sg){

      temp_res <- stats::cor.test(x_sg, y_sg,
                                  alternative = alternative,
                                  method = method,
                                  exact = exact,
                                  conf.level = conf_level,
                                  ...)

      tibble::tibble(
        !!names(temp_res$estimate) := temp_res$estimate,
        p.value = temp_res$p.value
      )
    }) %>% tibble::enframe(name = "y", value = "value") %>%
      tidyr::unnest(cols = value)


  }) %>% tibble::enframe(name = "x", value = "value") %>%
    tidyr::unnest(cols = value)

  cor_res

}

#' Proportional Hazards Regression Model analysis
#'
#' Fits a Cox proportional hazards regression model and extract log-rank test
#' results and tidy informations. log-rank test is just a special case of the
#' Cox model \code{(score test)}
#'
#' @param data a \code{data.frame} in which to interpret the variables named in
#'   the formula, or in the subset and the weights argument.
#' @param formula a formula object, with the response on the left of a ~
#'   operator, and the terms on the right. The response must be a survival
#'   object as returned by the Surv function.
#' @param exponentiate Logical indicating whether or not to exponentiate the the
#'   coefficient estimates. This is typical for logistic and multinomial
#'   regressions, but a bad idea if there is no log or logit link. Defaults:
#'   \code{TRUE}. Details see \code{\link[broom:tidy.coxph]{tidy}}
#' @param conf_level The confidence level to use for the confidence interval.
#'   Must be strictly greater than 0 and less than 1. Defaults to 0.95, which
#'   corresponds to a 95 percent confidence interval. Details see
#'   \code{\link[broom:tidy.coxph]{tidy}}. Default: \code{0.95}
#' @param ... Additional arguments passed to \code{\link[survival:coxph]{coxph}}
#' @return a list with two item. One for overall score test results and one for
#'   summarized information corresponding to each term in the RHS of formula.
#'   Details see \code{\link[broom:tidy.coxph]{tidy}}
#' @examples
#' survival_example_data <- readRDS(system.file(
#'   "extdata", "survival_example_data.rds", package = "yjkit"
#' ))
#' stat_cox_test( survival_example_data,
#'                survival::Surv(time, status) ~ ph.ecog + tt(age) )
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @references \itemize{\item
#'   \href{https://stats.stackexchange.com/questions/362381/logrank-p-value-for-2-groups}{logrank-p-value-for-2-groups}
#'    \cr \item
#'   \href{https://stats.stackexchange.com/questions/486806/the-logrank-test-statistic-is-equivalent-to-the-score-of-a-cox-regression-is-th}{the-logrank-test-statistic-is-equivalent-to-the-score-of-a-cox-regression-is-th}}
#' @export
stat_cox_test <- function(data, formula,
                          exponentiate = TRUE, conf_level = 0.95,
                          ...){

  stopifnot(inherits(data, "data.frame"))

  res <- survival::coxph(
    formula = formula,
    data = data, ...
  )

  list(
    sctest_res =  tibble::tibble(
      overall_var = deparse1(formula[[3]]),
      p.value = summary(res)$sctest["pvalue"]
    ),

    term_res = broom::tidy(res, exponentiate = exponentiate,
                           conf.int = TRUE,
                           conf.level = conf_level)
  )


}


# stat_between_test utility function --------------------------------------

stat_between_test_helper <- function(data, x, y, x_label = NULL, y_label = NULL,
                                     type = c("nonparametric", "parametric"),
                                     ...){

  stopifnot(inherits(data, "data.frame"))
  type <- match.arg(type)

  quo_x <- rlang::enquo(x)
  quo_y <- rlang::enquo(y)
  if (is.null(x_label)) x_label <- deparse1(rlang::enexpr(x))
  if (is.null(y_label)) y_label <- deparse1(rlang::enexpr(x))
  arg_dots <- rlang::enquos(...)

  test_data <- dplyr::mutate(
    data, ..x = !!quo_x, ..y = !!quo_y
  )

  x_unique_levels <- length(unique(stats::na.omit(test_data[["..x"]])))

  if(x_unique_levels < 2) stop("the unique values (omit missing values) for x in data should be at least 2")

  if (!is.numeric(test_data[["..y"]])) stop("y in data must be a numeric vector")

  if (type == "parametric") {

    if (x_unique_levels == 2) test_method <- "t.test"
    if (x_unique_levels > 2) test_method <- "aov"

  }

  if (type == "nonparametric") {

    if (x_unique_levels == 2) test_method <- "wilcox.test"
    if (x_unique_levels > 2) test_method <- "kruskal.test"

  }

  test_fn <- rlang::call2(
    "::", rlang::expr(stats), rlang::sym(test_method)
  )

  test_res <- rlang::expr(
    (!!test_fn)(..y ~ ..x, data = test_data, !!!arg_dots)
  )

  test_res <- broom::tidy( rlang::eval_tidy(test_res) )

  if (test_method == "aov") {
    test_res <- test_res %>%
      dplyr::filter(term == "..x") %>%
      dplyr::select(!term)
  }
  test_res <- test_res %>%
    dplyr::mutate(method = !!test_method,
                  term_y = y_label,
                  term_x = x_label) %>%
    dplyr::relocate(method, term_y, term_x, .before = 1)

  if (test_method %in% c("aov", "kruskal.test")) {
    names(test_res$statistic) <- switch(
      test_method,
      aov = "F",
      kruskal.test = "Chisq"
    )
  }
  test_res

}
