testthat::test_that(
  "test stat_between_test()",
  {
    testthat::expect_snapshot_value(
      stat_between_test(mtcars, factor(vs), mpg, type = "p"),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      stat_between_test(mtcars, factor(vs), mpg, type = "n"),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      stat_between_test(mtcars, factor(cyl), mpg, type = "p"),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      stat_between_test(mtcars, factor(cyl), mpg, type = "n"),
      style = "serialize"
    )
  }
)


testthat::test_that(
  "test stat_cor_test()",
  {
    testthat::expect_snapshot_value(
      stat_cor_test(mtcars),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      stat_cor_test(mtcars, cor_test = TRUE),
      style = "serialize"
    )

    # c("two.sided", "less", "greater") method = "spearman"
    testthat::expect_equal(
      stat_cor_test(mtcars, alternative = "two.sided"),
      stat_cor_test(mtcars, alternative = "two.sided", cor_test = TRUE)
    )

    testthat::expect_equal(
      stat_cor_test(mtcars, alternative = "less"),
      stat_cor_test(mtcars, alternative = "less", cor_test = TRUE)
    )

    testthat::expect_equal(
      stat_cor_test(mtcars, alternative = "greater"),
      stat_cor_test(mtcars, alternative = "greater", cor_test = TRUE)
    )

    # c("two.sided", "less", "greater") method = "pearson"
    testthat::expect_equal(
      stat_cor_test(mtcars, alternative = "two.sided", method = "pearson"),
      stat_cor_test(mtcars,
        alternative = "two.sided",
        method = "pearson", cor_test = TRUE
      )
    )

    testthat::expect_equal(
      stat_cor_test(mtcars, alternative = "less", method = "pearson"),
      stat_cor_test(mtcars,
        alternative = "less",
        method = "pearson", cor_test = TRUE
      )
    )

    testthat::expect_equal(
      stat_cor_test(mtcars, alternative = "greater", method = "pearson"),
      stat_cor_test(mtcars,
        alternative = "greater",
        method = "pearson", cor_test = TRUE
      )
    )
  }
)
