testthat::test_that(
  "test run_absolute()", {
    seg <- readRDS(system.file("extdata", "run_absolute_example_seg.rds",
                               package = "yjtools"))

    maf <- readRDS(system.file("extdata", "run_absolute_example_maf.rds",
                               package = "yjtools"))

    testthat::expect_output(
      run_absolute(
        seg = seg, maf = maf,
        results_dir = file.path(tempdir(),"results", "ABSOLUTE"),
        BPPARAM = BiocParallel::SerialParam()
      ),
      "ABSOLUTE algorithm Done.", fixed = TRUE
    )

    testthat::expect_snapshot_value(
      list.files(file.path(tempdir(),"results", "ABSOLUTE", "RunAbsolute")),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      list.files(file.path(tempdir(),"results", "ABSOLUTE",
                           "CreateReviewObject")),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      list.files(file.path(tempdir(),"results", "ABSOLUTE", "reviewed")),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      read.table(
        file.path(tempdir(),"results", "ABSOLUTE", "CreateReviewObject",
                  "SummarizeAbsolute.PP-calls_tab.txt"),
        sep = "\t", header = TRUE
      ),
      style = "serialize"
    )
    testthat::expect_snapshot_value(
      read.table(
        file.path(tempdir(),"results", "ABSOLUTE", "reviewed",
                  "ReviewAbsolute.YJ.ABSOLUTE.table.txt"),
        sep = "\t", header = TRUE
      ),
      style = "serialize"
    )


  }
)
