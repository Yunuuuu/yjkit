testthat::test_that(
  "test stat_between_test()", {

    testthat::expect_snapshot_value(
      stat_between_test(mtcars, vs, mpg, type = "p"),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      stat_between_test(mtcars, vs, mpg, type = "n"),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      stat_between_test(mtcars, cyl, mpg, type = "p"),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      stat_between_test(mtcars, cyl, mpg, type = "n"),
      style = "serialize"
    )

  }
)


testthat::test_that(
  "test stat_cor_test()", {

    testthat::expect_snapshot_value(
      stat_cor_test(mtcars),
      style = "serialize"
    )

  }
)


testthat::test_that(
  "test stat_cox_test()", {

    survival_example_data <- readRDS(
      system.file("extdata", "survival_example_data.rds", package = "yjtools")
    )
    testthat::expect_snapshot_value(
      stat_cox_test(survival_example_data,
                    survival::Surv(time, status) ~ ph.ecog + tt(age)),
      style = "serialize"
    )

  }
)
