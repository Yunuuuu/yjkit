#' Statistical Inference with Permutation
#'
#' @param data A data frame that can be coerced into a tibble.
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
#' @inheritDotParams infer::calculate -x -stat
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
#' @return A data.frame
#' @export
infer <- function(
    data, formula, stat, reps = 2000L, level = 0.95,
    direction = "two-sided", null = NULL, type = NULL,
    null_reps = reps, ci_reps = reps, ...,
    variables = NULL, response = NULL, explanatory = NULL, success = NULL,
    p = NULL, mu = NULL, med = NULL, sigma = NULL) {
    assert_pkg("infer")
    response <- rlang::enquo(response)
    explanatory <- rlang::enquo(explanatory)
    if (missing(formula)) {
        model <- infer::specify(data,
            response = !!response,
            explanatory = !!explanatory,
            success = success
        )
    } else {
        model <- infer::specify(data,
            formula = formula,
            response = !!response,
            explanatory = !!explanatory,
            success = success
        )
    }
    explanatory <- infer:::explanatory_name(model)
    if (length(explanatory) > 1L) {
    } else {
        stat <- match.arg(stat, infer:::implemented_stats)
    }
    # check if we need Declare a null hypothesis ---------------
    null_hypo <- quote(infer::hypothesize(
        model,
        null = null,
        p = p, mu = mu, med = med, sigma = sigma
    ))
    if (length(explanatory)) {
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
    # Calculating the observed statistic
    if (length(explanatory) > 1L) {
        obs_stat <- infer::fit(model)
        model <- eval(null_hypo)
        null_dist <- model
    } else {
        if (!any(stat == infer:::untheorized_stats)) {
            model <- eval(null_hypo)
        }
        obs_stat <- infer::calculate(x = model, stat = stat, ...)
        if (any(stat == infer:::untheorized_stats)) {
            null_dist <- eval(null_hypo)
        } else {
            null_dist <- model
        }
    }
    # generating the null distribution
    if (null_type == "permute") {
        variables <- rlang::enquo(variables)
        if (rlang::quo_is_null(variables)) {
            variables <- attr(model, "response")
        }
        null_dist <- infer::generate(null_dist,
            reps = null_reps, type = null_type,
            variables = !!variables
        )
    } else {
        null_dist <- infer::generate(null_dist,
            reps = null_reps, type = null_type
        )
    }
    if (length(explanatory) > 1L) {
        null_dist <- infer::fit(null_dist)
        # calculate Confidence intervals
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
