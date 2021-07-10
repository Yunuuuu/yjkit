#' Create a ggplot2 object with annotation for group or condition test results
#'
#' annotate group or condition comparisons in between-subjects designs in a
#' ggplot2 object.
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
#' @param format_label if \code{TRUE}, numeric label will be formatted by
#'   \code{\link{format_num}}. Default: \code{TRUE}
#' @param label_position the position where the label will added, can be a
#'   numeric vector of length two or a character vector with length one (one of
#'   \code{"lefttop", "righttop", "leftbottom", and "rightbottom"}) or length
#'   two (See below \code{label_justification}). Default: \code{c(0.02, 0.98)}.
#' @param label_justification label_justification must be a element-two (or one)
#'   character or numeric vector. Possible string values are: "left", "right",
#'   "centre", "center", "bottom", and "top". For numeric values, 0 means left
#'   (bottom) alignment and 1 means right (top) alignment. Default: \code{c(0,
#'   1)}
#' @param anno_statistic a scalar logical, Should test statistic be annotated in
#'   the label
#' @param label_sep see argument \code{label}. String to insert between the test
#'   statistic and test P-value.
#' @param anno_args  Other arguments passed on to \code{\link[ggplot2]{layer}}.
#'   These are often aesthetics, used to set an aesthetic to a fixed value, like
#'   colour = "red" or size = 3. They may also be parameters to the paired
#'   geom/stat. See \code{\link[ggplot2:geom_text]{geom_text}}
#' @param ... other arguments passed to statistical test function,
#'   \code{\link[stats]{wilcox.test}}, \code{\link[stats]{kruskal.test}},
#'   \code{\link[stats]{t.test}} and \code{\link[stats]{aov}}.
#' @return a ggplot2 object with annotation of group or condition comparisons
#'   test results
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#' gganno_between_test(mtcars, aes(factor(cyl), mpg)) +
#'   ggplot2::geom_boxplot()
#' @export
gganno_between_test <- function(data, mapping = NULL,
                                type = c("nonparametric", "parametric"),
                                format_label = TRUE,
                                label_position = c(0.02, 0.98),
                                label_justification = c(0, 1),
                                anno_statistic = TRUE,
                                label_sep = NULL,
                                anno_args = list(fontface = "bold"),
                                ...){

  test_aes <- c("x", "y") %in% names(mapping)
  if (!all(test_aes))  stop(
    paste0("lack of ", c("x", "y")[!test_aes], " aesthetic mappings")
  )

  type <- match.arg(type)

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

  ggplot2::ggplot(data = data, mapping = mapping)+
  rlang::exec(ggstat_anno_helper, geom = "text",
              label = label, format_label = format_label,
              label_position = label_position,
              label_justification = label_justification,
              label_sep = label_sep, !!!anno_args)
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
#' @param format_label if \code{TRUE}, numeric label will be formatted by
#'   \code{\link{format_num}}. Default: \code{TRUE}
#' @param label_position the position where the label will added, can be a
#'   numeric vector of length two or a character vector with length one (one of
#'   \code{"lefttop", "righttop", "leftbottom", and "rightbottom"}) or length
#'   two (See below \code{label_justification}). Default: \code{c(0.02, 0.98)}.
#' @param label_justification label_justification must be a element-two (or one)
#'   character or numeric vector. Possible string values are: "left", "right",
#'   "centre", "center", "bottom", and "top". For numeric values, 0 means left
#'   (bottom) alignment and 1 means right (top) alignment. Default: \code{c(0,
#'   1)}.
#' @param anno_statistic a scalar logical, Should test statistic added in label
#' @param label_sep see argument \code{label}. String to insert between the test
#'   statistic and test P-value.
#' @param anno_args  Other arguments passed on to \code{\link[ggplot2]{layer}}.
#'   These are often aesthetics, used to set an aesthetic to a fixed value, like
#'   colour = "red" or size = 3. They may also be parameters to the paired
#'   geom/stat. See \code{\link[ggplot2:geom_text]{geom_text}}
#' @param ... Other parameters for \code{\link[stats]{cor.test}}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#'   gganno_cor_test(mtcars, aes(disp, mpg)) +
#'     ggplot2::geom_point()+
#'     ggplot2::geom_smooth()
#' @export
gganno_cor_test <- function(data, mapping = NULL,
                            use = c("pairwise.complete.obs",
                                    "everything", "all.obs", "complete.obs",
                                    "na.or.complete"),
                            method = c("spearman", "pearson", "kendall"),
                            alternative = c("two.sided", "less", "greater"),
                            exact = TRUE, conf_level = 0.95,
                            format_label = TRUE,
                            label_position = c(0.02, 0.98),
                            label_justification = c(0, 1),
                            anno_statistic = TRUE,
                            label_sep = NULL,
                            anno_args = list(fontface = "bold"),
                            ...){

  test_aes <- c("x", "y") %in% names(mapping)
  if (!all(test_aes))  stop(
    paste0("lack of ", c("x", "y")[!test_aes], " aesthetic mappings")
  )

  use <- match.arg(use)
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

  ggplot2::ggplot(data = data, mapping = mapping)+
    rlang::exec(ggstat_anno_helper, geom = "text",
                label = label, format_label = format_label,
                label_position = label_position,
                label_justification = label_justification,
                label_sep = label_sep, !!!anno_args)

}


# ggstat_anno_helper ------------------------------------------------

#' Create an annotation layer with text annotation
#'
#' This function Creates an ggplot2 object with annotation.
#'
#' @param label the character which will be added to the plot, will be connected
#'   by \code{paste(names(label), label, collapse = label_sep, sep = ": ")}
#' @param format_label if \code{TRUE}, numeric label will be formatted by
#'   \code{\link{format_num}}. Default: \code{TRUE}
#' @param label_position the position where the label will added, can be a
#'   numeric vector of length two or a character vector with length one (one of
#'   \code{"lefttop", "righttop", "leftbottom", and "rightbottom"}) or length
#'   two (See below \code{label_justification}). Default: \code{c(0.02, 0.98)}.
#' @param label_justification label_justification must be a element-two (or one)
#'   character or numeric vector. Possible string values are: "left", "right",
#'   "centre", "center", "bottom", and "top". For numeric values, 0 means left
#'   (bottom) alignment and 1 means right (top) alignment. Default: \code{c(0,
#'   1)}
#' @param label_sep see argument \code{label}. String to insert between each of
#'   the \code{label}.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}.
#'   These are often aesthetics, used to set an aesthetic to a fixed value, like
#'   colour = "red" or size = 3. They may also be parameters to the paired
#'   geom/stat. See \code{\link[ggplot2:geom_text]{geom_text}}
#' @return a ggplot2 annotation layer
# @examples
#   ggstat_anno(label = "a", label_position = c(0.5, 0.5))
#' @author Yun \email{yunyunpp96@@outlook.com}
ggstat_anno_helper <- function(label, format_label = TRUE,
                               label_position = c(0.02, 0.98),
                               label_justification = c(0, 1),
                               label_sep = NULL,
                               ...){

  if(identical(length(label_position), 1L) && is.null(label_sep)) {

    label_sep <- "; "

  }

  if(identical(length(label_position), 2L) && is.null(label_sep)) {

    label_sep <- "\n"

  }

  if (format_label && is.numeric(label)) {
    label_format <- format_num(label)
  } else {
    label_format <- as.character(label)
  }

  label <- stringr::str_c(names(label), label_format, sep = ": ",
                          collapse = label_sep)

  ggannotate_npc(geom = "text",
                 label = label,
                 label_position = label_position,
                 label_justification,
                 ...)

}



