#' group or condition comparisons in between-subjects designs.
#'
#' Create a ggplot2 object with annotation for group or condition test results
#'
#' @param data Default dataset to use for plot. If not specified, must be
#'   supplied in each layer added to the plot.
#' @param mapping Default list of aesthetic mappings to use for plot. x and y
#'   aesthetic must in.
#' @param type if \code{type = "nonparametric"},
#'   \code{\link[stats]{wilcox.test}} or \code{\link[stats]{kruskal.test}} is
#'   used. \cr if \code{type = "parametric"}, \code{\link[stats]{t.test}} or
#'   \code{\link[stats]{aov}} is used.\cr which method is used depends on the
#'   unique values of \code{x} in \code{data}
#' @param format_label if \code{TRUE}, will be formatted by
#'   \code{\link{format_num}}. Default: \code{TRUE}
#' @param label_position the position where the label will added, can be a
#'   numeric vector of length 2 or a character vector \code{(one of "title",
#'   "subtitle", "caption", and "tag")} of length one.
#' @param label_justification These can either be a number between 0
#'   (\code{right/bottom}) and 1 (\code{top/left}) or a character
#'   ("\code{left}", "\code{middle}", "\code{right}", "\code{bottom}",
#'   "\code{center}", "\code{top}"). There are two special alignments:
#'   "\code{inward}" and "\code{outward}". Inward always aligns text towards the
#'   center, and outward aligns it away from the center.
#' @param anno_statistic a scalar logical, Should test statistic added in label
#' @param label_sep see argument \code{label}. String to insert between the test
#'   statistic and test P-value.
#' @param ggtheme a ggplot2 theme object. Provided as a ggplot2 theme object or
#'   a function (can be a function name as a string) implemented as a ggplot2
#'   theme object.
#' @param anno_args additional arguments control the \code{label}. See
#'   \code{\link[grid]{textGrob}} and \code{\link[ggplot2]{element_text}}
#' @param ... other arguments passed to statistical test function,
#'   \code{\link[stats]{wilcox.test}}, \code{\link[stats]{kruskal.test}},
#'   \code{\link[stats]{t.test}} and \code{\link[stats]{aov}}
#' @return a ggplot2 object with annotation of group or condition comparisons
#'   test results
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#'   gganno_between_test(mtcars, aes(cyl, mpg)) +
#'   geom_boxplot()
#' @export
gganno_between_test <- function(data, mapping = aes(x, y, ...),
                                type = c("nonparametric", "parametric"),
                                format_label = TRUE,
                                label_position = c(0.02, 1),
                                label_justification = NULL,
                                anno_statistic = TRUE,
                                label_sep = NULL,
                                ggtheme = ggthemes::theme_few(),
                                anno_args = list(fontface = "bold"),
                                ...){

  test_aes <- c("x", "y") %in% names(mapping)
  if (!all(test_aes))  stop(
    paste0("lack of ", c("x", "y")[!test_aes], " aesthetic mappings")
  )

  test_res <- stat_between_test_helper(data = data,
                                       x = !!mapping$x,
                                       y = !!mapping$y,
                                       type = type,
                                       ...)


  if (anno_statistic) {
    label <- c(test_res$statistic, test_res$p.value) %>%
      rlang::set_names(
        nm = c(names(test_res$statistic), "P-value")
      )
  } else {
    label <- c(`P-value` = test_res$p.value)
  }

  rlang::exec(gganno_text, data = data, mapping = mapping,
              label = label, format_label = format_label,
              label_position = label_position,
              label_justification = label_justification,
              label_sep = label_sep, ggtheme, !!!anno_args)
}


#' annotate Correlation coefficient in a ggplot2 object
#'
#' Create a ggplot2 object with annotation for Correlation coefficient
#'
#' @param data Default dataset to use for plot. If not specified, must be
#'   supplied in each layer added to the plot.
#' @param mapping Default list of aesthetic mappings to use for plot. x and y
#'   aesthetic must in.
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
#' @param format_label if \code{TRUE}, will be formatted by
#'   \code{\link{format_num}}. Default: \code{TRUE}
#' @param label_position the position where the label will added, can be a
#'   numeric vector of length 2 or a character vector \code{(one of "title",
#'   "subtitle", "caption", and "tag")} of length one.
#' @param label_justification These can either be a number between 0
#'   (\code{right/bottom}) and 1 (\code{top/left}) or a character
#'   ("\code{left}", "\code{middle}", "\code{right}", "\code{bottom}",
#'   "\code{center}", "\code{top}"). There are two special alignments:
#'   "\code{inward}" and "\code{outward}". Inward always aligns text towards the
#'   center, and outward aligns it away from the center.
#' @param anno_statistic a scalar logical, Should test statistic added in label
#' @param label_sep see argument \code{label}. String to insert between the test
#'   statistic and test P-value.
#' @param ggtheme a ggplot2 theme object. Provided as a ggplot2 theme object or
#'   a function (can be a function name as a string) implemented as a ggplot2
#'   theme object.
#' @param anno_args additional arguments control the \code{label}. See
#'   \code{\link[grid]{textGrob}} and \code{\link[ggplot2]{element_text}}
#' @param ... Other parameters for \code{\link[stats]{cor.test}}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#'   gganno_cor_test(mtcars, aes(disp, mpg)) +
#'     ggplot2::geom_point()+
#'     ggplot2::geom_smooth()
#' @export
gganno_cor_test <- function(data, mapping = aes(x, y, ...),
                            use = c("pairwise.complete.obs",
                                    "everything", "all.obs", "complete.obs",
                                    "na.or.complete"),
                            method = c("spearman", "pearson", "kendall"),
                            alternative = c("two.sided", "less", "greater"),
                            exact = TRUE, conf_level = 0.95,
                            format_label = TRUE,
                            label_position = c(0.02, 1),
                            label_justification = NULL,
                            anno_statistic = TRUE,
                            label_sep = NULL,
                            ggtheme = ggthemes::theme_few(),
                            anno_args = list(fontface = "bold"),
                            ...){

  test_aes <- c("x", "y") %in% names(mapping)
  if (!all(test_aes))  stop(
    paste0("lack of ", c("x", "y")[!test_aes], " aesthetic mappings")
  )

  use = match.arg(use)
  method <- match.arg(method)
  alternative <- match.arg(alternative)

  test_res <- rlang::expr(
    stats::cor.test(formula = ~ !!(mapping$x) + !!(mapping$y),
                    data = data,
                    use = use,
                    method = method,
                    alternative = alternative,
                    exact = exact,
                    conf.level = conf_level, ...)
  )
  test_res <- rlang::eval_tidy(test_res)


  if (anno_statistic) {
    label <- c(test_res$estimate, test_res$p.value) %>%
      rlang::set_names(
        nm = c(names(test_res$estimate), "P-value")
      )
  } else {
    label <- c(`P-value` = test_res$p.value)
  }

  rlang::exec(gganno_text, data = data, mapping = mapping,
              label = label, format_label = format_label,
              label_position = label_position,
              label_justification = label_justification,
              label_sep = label_sep, ggtheme, !!!anno_args)
}

