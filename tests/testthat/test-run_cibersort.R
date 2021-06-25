testthat::test_that(
  "test run_cibersort()", {

    testthat::expect_snapshot_value(
      run_cibersort(
        readRDS(
          system.file("extdata", "run_cibersort_example_mixture_data.rds",
                      package = "yjkit")
        ),
        quantile_norm = TRUE, perm = 50,
        absolute = FALSE,
        abs_method = "sig.score"
      ),
      style = "serialize"
    )
  }

)
