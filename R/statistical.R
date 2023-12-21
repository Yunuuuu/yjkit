#' Between group test
#'
#' Function to calculate statistics and P-value for both parametric and
#' nonparametric tests.
#'
#' @param data a data.frame object containing the variables in the x and y or
#'   NULL.
#' @param x,y if \code{data} supplied, \code{x} and \code{y} are quoted
#'   arguments, can be a list, if \code{x} or \code{y} is a list object, each
#'   element in \code{x} and in \code{y} will be paired to implement test; if
#'   \code{data} is NULL, \code{x} and \code{y} should be a vector or a list
#'   with atmomic vector in the same length.
#' @param type if \code{type = "nonparametric"},
#'   \code{\link[stats]{wilcox.test}} or \code{\link[stats]{kruskal.test}} is
#'   used. \cr if \code{type = "parametric"}, \code{\link[stats]{t.test}} or
#'   \code{\link[stats]{aov}} is used.\cr which method is used depends on the
#'   unique values of \code{x} in \code{data}
#' @param ... other arguments passed to statistical test function,
#'   \code{\link[stats]{wilcox.test}}, \code{\link[stats]{kruskal.test}},
#'   \code{\link[stats]{t.test}} and \code{\link[stats]{aov}}
#' @section Quasiquotation: if \code{data} supplied, \code{x} and \code{y} in
#'   \code{stat_between_test()} are \code{\link[rlang:nse-defuse]{quoted
#'   arguments}} This means that its inputs are quoted to be evaluated in the
#'   context of the data. This makes it easy to work with variables from the
#'   data frame because you can name those directly. The flip side is that you
#'   have to use \code{\link[rlang:nse-force]{quasiquotation}} to program with
#'   \code{stat_between_test()}. See a tidy evaluation tutorial such as the
#'   \href{https://dplyr.tidyverse.org/articles/programming.html}{dplyr
#'   programming vignette} to learn more about these techniques.
#'
#' @return a tibble of the test results where each row match each pair of x and
#'   y. (see \code{\link[generics:tidy]{tidy}})
#' @details if \code{type = "nonparametric"}, \code{\link[stats]{wilcox.test}}
#'   or \code{\link[stats]{kruskal.test}} is used. \cr if \code{type =
#'   "parametric"}, \code{\link[stats]{t.test}} or \code{\link[stats]{aov}} is
#'   used.\cr The final choice of which method is used depends on the unique
#'   values of \code{x} term.
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#' stat_between_test(mtcars, factor(vs), mpg, type = "p")
#' stat_between_test(mtcars, factor(vs), mpg, type = "n")
#' stat_between_test(mtcars, factor(cyl), mpg, type = "p")
#' stat_between_test(mtcars, factor(cyl), mpg, type = "n")
#' stat_between_test(
#'     x = lapply(mtcars[c(2, 8)], as.factor),
#'     y = mtcars[c(1, 3, 4)], type = "n"
#' )
#' stat_between_test(
#'     x = factor(mtcars[["vs"]]),
#'     y = mtcars["mpg"], type = "n"
#' )
#' @export
stat_between_test <- function(data = NULL, x, y,
                              type = c("nonparametric", "parametric"),
                              ...) {
    type <- match.arg(type)

    if (!is.null(data)) {
        stopifnot(inherits(data, "data.frame"))
        quo_list_x <- rlang::enquos(x)
        label_list_x <- rlang::enexprs(x)
        quo_list_y <- rlang::enquos(y)
        label_list_y <- rlang::enexprs(y)

        y_list <- lapply(seq_along(quo_list_y), function(y_i) {
            x_list <- lapply(seq_along(quo_list_x), function(x_i) {
                stat_between_test_helper(
                    data = data,
                    x = !!quo_list_x[[x_i]],
                    y = !!quo_list_y[[y_i]],
                    x_label = deparse1(label_list_x[[x_i]]),
                    y_label = deparse1(label_list_y[[y_i]]),
                    type = type, ...
                )
            })

            dplyr::bind_rows(x_list)
        })

        res <- dplyr::bind_rows(y_list)
    } else {
        label_list_x <- rlang::enexprs(x)
        label_list_y <- rlang::enexprs(y)

        # check the right type of x
        if (is.character(x) || is.factor(x)) {
            x <- list(x = x)
            label_list_x <- lapply(label_list_x, deparse1)
        } else if (is.list(x)) {
            if (!all(vapply(x, function(x) is.character(x) || is.factor(x), logical(1)))) {
                stop("x should be a character or a factor or a list ",
                    "of character or factor vector with same length.",
                    call. = FALSE
                )
            }

            label_list_x <- names(x)
        } else {
            stop("unsupported type of ", typeof(x), " x, ",
                "x should be a character or a factor or a list ",
                "of character or factor vector with same length.",
                call. = FALSE
            )
        }


        # check the right type of y
        if (is.numeric(y)) {
            y <- list(y = y)
            label_list_y <- lapply(label_list_y, deparse1)
        } else if (is.list(y)) {
            if (!all(vapply(y, is.numeric, logical(1)))) {
                stop("y should be a numeric vector or a list ",
                    "of numeric vector with same length.",
                    call. = FALSE
                )
            }

            label_list_y <- names(y)
        } else {
            stop("unsupported type of ", typeof(y), " y, ",
                "y should be a numeric vector or a list ",
                "of numeric vector with same length.",
                call. = FALSE
            )
        }

        # check the same length of x and y
        if (!identical(
            length(unique(vapply(c(x, y), length, integer(1)))), 1L
        )) {
            stop("the length of each vector in x and in y should be the same",
                call. = FALSE
            )
        }

        y_list <- lapply(seq_along(y), function(y_i) {
            x_list <- lapply(seq_along(x), function(x_i) {
                stat_between_test_helper(
                    data = NULL,
                    x = x[[x_i]],
                    y = y[[y_i]],
                    x_label = label_list_x[[x_i]],
                    y_label = label_list_y[[y_i]],
                    type = type, ...
                )
            })
            dplyr::bind_rows(x_list)
        })

        res <- dplyr::bind_rows(y_list)
    }

    res
}

#' Correlation analysis
#'
#' Function to calculate Correlation coefficient and P-value between the columns
#' of x and the columns of y, using one of Pearson's product moment correlation
#' coefficient, Kendall's tau or Spearman's rho.
#'
#' @param x,y data.frame object. x and y must have the same \code{nrow}, every
#'   column in x or y will be coerced to a numeric vector.
#' @param use an optional character string giving a method for computing
#'   covariances in the presence of missing values. This must be (an
#'   abbreviation of) one of the strings \code{"everything"}, \code{"all.obs"},
#'   \code{"complete.obs"}, \code{"na.or.complete"}, or
#'   \code{"pairwise.complete.obs"}. Default: \code{"pairwise.complete.obs"}.
#'   See \code{\link[stats]{cor}}.
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
#' @param cor_test a logical value indicates whether the estimation of P-value
#'   should use \code{\link[stats]{cor.test}} instead. If FALSE, the estimation
#'   of P-value will use [stats::pt()] (for spearman or pearson) or
#'   [stats::pnorm()] (for kendall). Default: \code{FALSE}.
#' @param padj_method See [p.adjust][stats::p.adjust()] `method`.
#' @param ... Other parameters passed to \code{\link[stats]{cor.test}}.
#' @return a tibble of the Correlation results. (see
#'   \code{\link[stats]{cor.test}}).
#' @details See \code{\link[stats]{cor.test}} and \code{\link[stats]{cor}}.
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @examples
#' stat_cor_test(mtcars)
#' stat_cor_test(mtcars, cor_test = TRUE)
#' @export
stat_cor_test <- function(x, y = NULL,
                          use = c(
                              "pairwise.complete.obs",
                              "everything", "all.obs", "complete.obs",
                              "na.or.complete"
                          ),
                          method = c("spearman", "pearson", "kendall"),
                          alternative = c("two.sided", "less", "greater"),
                          exact = TRUE, cor_test = FALSE,
                          padj_method = "fdr",
                          ...) {
    assert_data_frame(x)
    assert_bool(cor_test)

    if (is.null(y)) {
        y <- x
    } else {
        stopifnot(inherits(y, "data.frame"))
        if (!identical(nrow(x), nrow(y))) {
            stop(
                "x and y must have same number of rows",
                call. = FALSE
            )
        }
    }

    use <- match.arg(use)
    method <- match.arg(method)
    alternative <- match.arg(alternative)

    x <- dplyr::mutate(x, dplyr::across(.fns = as.double))
    y <- dplyr::mutate(y, dplyr::across(.fns = as.double))

    cor_name <- switch(method,
        pearson = "cor",
        kendall = "tau",
        spearman = "rho"
    )

    if (cor_test) { # too slowly for large data

        res <- purrr::imap_dfr(x, function(data_x, name_x) {
            each_x_cor_y <- purrr::imap_dfr(y, function(data_y, name_y) {
                temp_res <- stats::cor.test(data_x, data_y,
                    alternative = alternative,
                    method = method,
                    exact = exact,
                    ...
                )

                tibble::tibble(
                    y = name_y,
                    !!cor_name := temp_res$estimate,
                    pvalue = temp_res$p.value
                )
            })

            dplyr::mutate(each_x_cor_y, x = !!name_x, .before = 1)
        })

        res[[cor_name]] <- unname(res[[cor_name]])
    } else {
        cor_res <- stats::cor(x, y, use = use, method = method) %>%
            tibble::as_tibble(rownames = ".x.", .name_repair = "minimal") %>%
            tidyr::pivot_longer(
                cols = !dplyr::all_of(".x."),
                names_to = ".y.",
                values_to = "cor"
            ) %>%
            rlang::set_names(nm = c("x", "y", cor_name))

        if (identical(use, "pairwise.complete.obs")) {
            x_matrix <- as.matrix(x)
            y_matrix <- as.matrix(y)

            n_matrix <- t(!is.na(x_matrix)) %*% (!is.na(y_matrix))
            n <- n_matrix[as.matrix(cor_res[c("x", "y")])]
        } else if (identical(use, "complete.obs") || identical(use, "na.or.complete")) {
            bind_x_and_y <- tidyr::drop_na(
                dplyr::bind_cols(x, y, .name_repair = "minimal")
            )
            n <- nrow(bind_x_and_y)
        } else {
            n <- nrow(x)
        }

        res <- dplyr::mutate(
            cor_res,
            pvalue = cor_to_p(
                .data[[cor_name]],
                n = !!n,
                method = !!method,
                alternative = !!alternative
            )
        )
    }

    dplyr::mutate(
        res,
        padjust = stats::p.adjust(.data[["pvalue"]], method = !!padj_method)
    )
}

# transform correlation coefficient to P-value
cor_to_p <- function(cor, n, method, alternative) {
    method <- match.arg(method, c("spearman", "pearson", "kendall"))
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

    if (identical(method, "kendall")) {
        warning("Estimation of P-value for Kendall's correlation ",
            "is not perfectly correct. Please try cor_test = TRUE.",
            call. = FALSE
        )
        statistic <- (3 * cor * sqrt(n * (n - 1))) / sqrt(2 * (2 * n + 5))
        p_fn <- "pnorm"
    } else {
        statistic <- cor * sqrt((n - 2) / (1 - cor^2))
        p_fn <- "pt"
    }

    if (!identical(alternative, "greater")) {
        lower.tail <- TRUE
    } else {
        lower.tail <- FALSE
    }

    if (identical(alternative, "two.sided")) statistic <- -abs(statistic)

    p_expr <- rlang::call2(p_fn,
        q = statistic,
        lower.tail = lower.tail,
        .ns = "stats"
    )

    if (!identical(method, "kendall")) {
        p_expr <- rlang::call_modify(p_expr, df = n - 2)
    }

    p <- rlang::eval_tidy(p_expr)

    if (identical(alternative, "two.sided")) p <- 2 * p

    p
}

# stat_between_test utility function --------------------------------------
stat_between_test_helper <- function(data = NULL, x, y,
                                     x_label = NULL, y_label = NULL,
                                     type = c("nonparametric", "parametric"),
                                     ...) {
    if (is.null(x_label)) x_label <- deparse1(rlang::enexpr(x))
    if (is.null(y_label)) y_label <- deparse1(rlang::enexpr(x))

    if (!is.null(data)) {
        stopifnot(inherits(data, "data.frame"))
        quo_x <- rlang::enquo(x)
        quo_y <- rlang::enquo(y)
        test_data <- dplyr::mutate(
            data,
            ..x = !!quo_x, ..y = !!quo_y
        )

        if (!(is.character(test_data[["..x"]]) || is.factor(test_data[["..x"]]))) {
            stop("x in data should be a type of character or factor but not a ",
                "type of ", typeof(test_data[["..x"]]),
                call. = FALSE
            )
        }

        x_unique_levels <- length(unique(stats::na.omit(test_data[["..x"]])))

        if (x_unique_levels < 2) stop("the unique values (omit missing values) for x in data should be at least 2", call. = FALSE)

        if (!is.numeric(test_data[["..y"]])) stop("y in data must be a numeric vector", call. = FALSE)
    } else {
        if (!identical(length(x), length(y))) {
            stop("the length of vector x and y should be the same",
                call. = FALSE
            )
        }

        if (!(is.character(x) || is.factor(x))) {
            stop("x should be a type of character or factor but not a ",
                "type of ", typeof(x),
                call. = FALSE
            )
        }

        x_unique_levels <- length(unique(stats::na.omit(x)))
        if (x_unique_levels < 2) stop("the unique values (omit missing values) for x should be at least 2", call. = FALSE)

        if (!is.numeric(y)) {
            stop("y should be a numeric vector but not a type of ",
                typeof(y),
                call. = FALSE
            )
        }

        test_data <- tibble::tibble(..x = x, ..y = y)
    }

    type <- match.arg(type)
    arg_dots <- rlang::enquos(...)

    if (type == "parametric") {
        if (x_unique_levels == 2) test_method <- "t.test"
        if (x_unique_levels > 2) test_method <- "aov"
    }

    if (type == "nonparametric") {
        if (x_unique_levels == 2) test_method <- "wilcox.test"
        if (x_unique_levels > 2) test_method <- "kruskal.test"
    }

    test_expr <- rlang::call2(
        test_method,
        ..y ~ ..x,
        data = test_data, !!!arg_dots,
        .ns = "stats"
    )

    test_res <- broom::tidy(rlang::eval_tidy(test_expr))

    if (test_method == "aov") {
        test_res <- test_res %>%
            dplyr::filter(.data$term == "..x") %>%
            dplyr::select(!dplyr::all_of("term"))
    }
    test_res <- test_res %>%
        dplyr::mutate(
            method = !!test_method,
            term_y = !!y_label,
            term_x = !!x_label
        ) %>%
        dplyr::relocate(dplyr::all_of("method"),
            dplyr::all_of("term_y"),
            dplyr::all_of("term_x"),
            .before = 1
        )

    if (test_method %in% c("aov", "kruskal.test")) {
        names(test_res$statistic) <- switch(test_method,
            aov = "F",
            kruskal.test = "Chisq"
        )
    }
    test_res
}
