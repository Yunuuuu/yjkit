#' Statistical Inference with Permutation
#'
#' @param data A data frame that can be coerced into a tibble.
#' @param paired A column in the data to specify the observation (sample id), if
#' not `NULL`, paired comparison will be applied.
#' @param variables A set of unquoted column names in the data to permute
#' (independently of each other). Defaults to only the response variable. Note
#' that any derived effects that depend on these columns (e.g., interaction
#' effects) will also be affected.
#' @param reps Fast way to set both `null_reps` and `ci_reps`.
#' @param null_reps Number of times to calculate null distribution.
#' @param ci_reps Number of bootstrap times to calculate confidence interval.
#' @inheritParams infer::get_ci
#' @inheritParams infer::calculate
#' @inheritParams infer::hypothesize
#' @inheritParams infer::get_pvalue
#' @inheritParams infer::specify
#' @inheritParams infer::generate
#' @param ... Other arguments passed to [infer::calculate].
#'  - `order`: A string vector of specifying the order in which the levels of
#'    the explanatory variable should be ordered for subtraction (or division
#'    for ratio-based statistics), where `order = c("first", "second")` means
#'    `("first" - "second")`, or the analogue for ratios. Needed for inference
#'    on difference in means, medians, proportions, ratios, t, and z statistics.
#'  - `...`: To pass options like `na.rm = TRUE` into functions like [mean()],
#'    [sd()], etc. Can also be used to supply hypothesized null values for the
#'    "t" statistic or additional arguments to [stats::chisq.test()].
#' @seealso
#' <https://infer.tidymodels.org/articles/observed_stat_examples.html>
#' @examples
#' data(gss, package = "infer")
#' infer(gss, response = hours, stat = "mean", mu = 40)
#' infer(gss, response = hours, stat = "t", mu = 40)
#' infer(gss, response = hours, stat = "median", med = 40)
#' infer(gss, response = sex, success = "female", stat = "prop", p = .5)
#' infer(gss, response = sex, success = "female", stat = "z", p = .5)
#' infer(gss, college ~ sex,
#'     success = "no degree",
#'     stat = "diff in props",
#'     order = c("female", "male")
#' )
#' infer(gss, hours ~ age + college, variables = c(age, college))
#' gss$hours_previous <- gss$hours + 5 - rpois(nrow(gss), 4.8)
#' gss$.id <- seq_len(nrow(gss))
#' gss_paired <- tidyr::pivot_longer(gss, cols = c(hours, hours_previous))
#' infer(gss_paired, value ~ name,
#'     stat = "mean", paired = .id
#' )
#' @return A [data.frame]
#' @export
infer <- function(
    data, formula, stat, paired = NULL, reps = 2000L, level = 0.95,
    direction = "two-sided", null = NULL, type = NULL,
    null_reps = reps, ci_reps = reps, ...,
    variables = NULL, response = NULL, explanatory = NULL, success = NULL,
    p = NULL, mu = NULL, med = NULL, sigma = NULL) {
    assert_pkg("infer")
    assert_data_frame(data)
    response <- rlang::enexpr(response)
    explanatory <- rlang::enexpr(explanatory)
    formula <- parse_formula(
        rlang::maybe_missing(formula),
        response, explanatory
    )
    paired <- rlang::enquo(paired)
    assert_paired_stat(paired, stat)

    # Specify response and explanatory variables --------------
    model <- specify_infer(data,
        paired = paired,
        response = formula$response,
        explanatory = formula$explanatory,
        success = success
    )

    # prepare arguments ---------------
    explanatory <- infer:::explanatory_name(model)
    if (length(explanatory) > 1L) {
        # always use infer::fit to calculate observed statistical
    } else {
        stat <- infer:::check_calculate_stat(stat)
    }
    if (!rlang::quo_is_null(paired)) {
        null <- null %||% "paired independence"
        null_type <- "permute"
    } else if (length(explanatory)) {
        null <- null %||% "independence"
        null_type <- "permute"
    } else {
        null <- null %||% "point"
        if (identical(attr(model, "response_type"), "factor")) {
            null_type <- "draw"
        } else {
            null_type <- "bootstrap"
        }
    }

    # prepare null hypothesis ---------------
    null_model <- infer::hypothesize(
        model,
        null = null,
        p = p, mu = mu, med = med, sigma = sigma
    )

    if (length(explanatory) > 1L) {
        # Calculating the observed statistic ---------
        obs_stat <- infer::fit(model)
    } else {
        # for theorized stats, we should declare null hypothesis first ----
        if (!any(stat == infer:::untheorized_stats)) {
            model <- null_model
        }
        # Calculating the observed statistic -----------------
        obs_stat <- infer::calculate(x = model, stat = stat, ...)
    }

    # generating the null distribution
    if (null_type == "permute") {
        variables <- rlang::enquo(variables)
        if (rlang::quo_is_null(variables)) {
            variables <- attr(null_model, "response")
        }
        null_dist <- infer::generate(null_model,
            reps = null_reps, type = null_type,
            variables = !!variables
        )
    } else {
        null_dist <- infer::generate(null_model,
            reps = null_reps, type = null_type
        )
    }
    # calculate Confidence intervals
    if (length(explanatory) > 1L) {
        null_dist <- infer::fit(null_dist)
        ci <- infer::get_ci(
            x = null_dist,
            level = level, type = type,
            point_estimate = obs_stat
        )
    } else {
        null_dist <- infer::calculate(x = null_dist, stat = stat, ...)
        boot_dist <- infer::generate(model,
            reps = ci_reps, type = attr(model, "type")
        )
        boot_dist <- infer::calculate(x = boot_dist, stat = stat, ...)
        # for small samples, boot_dist may contain `NA` values
        # waiting for `infer` upstream to fix
        # Now, we just remove the NA values before calculate CI for termporary
        # fix
        boot_dist <- boot_dist[!is.na(boot_dist$stat), ]
        # calculate Confidence intervals
        ci <- infer::get_ci(x = boot_dist, level = level, type = type)
    }
    # calculate P-value
    pvalue <- infer::get_pvalue(null_dist, obs_stat, direction = direction)
    if (length(explanatory) > 1L) {
        Reduce(function(x, y) {
            merge(x, y, by = "term", all = TRUE)
        }, list(obs_stat, pvalue, ci))
    } else {
        cbind(obs_stat, pvalue, ci)
    }
}

specify_infer <- function(data, paired, response, explanatory, ...) {
    if (!rlang::quo_is_null(paired)) {
        # parse paired formula ------------------------------------
        if (is.null(response)) {
            cli::cli_abort(
                "Please supply a response variable that is not {.code NULL}."
            )
        }
        if (n_vars(response) > 1L) {
            cli::cli_abort(
                "Only one response variable can be supplied for paired comparison"
            )
        }
        if (is.null(explanatory)) {
            cli::cli_abort(
                "Please supply a explanatory variable that is not {.code NULL} for paired comparison"
            )
        }
        if (n_vars(explanatory) > 1L) {
            cli::cli_abort(
                "Only one explanatory variable can be supplied for paired comparison"
            )
        }
        groups <- factor(pull(data, var = !!explanatory))
        n_group <- nlevels(groups)
        if (n_group != 2L) {
            cli::cli_abort(c(
                "Paired comparison can only be applied for two groups",
                i = "You have supplied {n_group} group{?s}"
            ))
        }
        data <- tidyr::pivot_wider(data,
            id_cols = !!paired,
            names_from = !!explanatory,
            values_from = !!response
        )
        group1 <- levels(groups)[1L]
        group2 <- levels(groups)[2L]
        data$.diff <- data[[group2]] - data[[group1]]
        response <- quote(.diff)
        infer::specify(x = data, response = !!response, ...)
    } else {
        infer::specify(
            x = data,
            response = !!response,
            explanatory = !!explanatory, ...
        )
    }
}

parse_formula <- function(formula, response, explanatory, call = rlang::caller_env()) {
    if (!rlang::is_missing(formula)) {
        if (!rlang::is_formula(formula)) {
            cli::cli_abort("{.arg formula} must be a {.cls formula}",
                call = call
            )
        }
        lhs <- rlang::f_lhs(formula)
        rhs <- rlang::f_rhs(formula)
        if (!is.null(lhs)) {
            response <- lhs
        }
        if (!is.null(rhs)) {
            explanatory <- rhs
        }
    }
    list(response = response, explanatory = explanatory)
}

n_vars <- function(x) {
    length(all.vars(x))
}

assert_paired_stat <- function(paired, stat, call = rlang::caller_env()) {
    if (!rlang::quo_is_null(paired)) {
        paired_stats <- c("mean", "median", "sum", "sd")
        if (!any(stat == paired_stats)) {
            cli::cli_abort(
                "Paired comparison can only applied for stat: {.val {paired_stats}}",
                call = call
            )
        }
    }
}
