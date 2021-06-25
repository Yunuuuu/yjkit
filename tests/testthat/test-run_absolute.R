testthat::test_that(
  "test run_absolute()", {
    seg <- readRDS(system.file("extdata", "run_absolute_example_seg.rds",
                               package = "yjkit"))

    maf <- readRDS(system.file("extdata", "run_absolute_example_maf.rds",
                               package = "yjkit"))

    testthat::expect_output(
      run_absolute(
        seg = seg, maf = maf,
        results_dir = file.path(tempdir(),"results", "ABSOLUTE"),
        BPPARAM = BiocParallel::SnowParam(
          BiocParallel::snowWorkers()
        )
      ),
      "-> Done.", fixed = TRUE
    )

    testthat::expect_snapshot_value(
      list.files(file.path(tempdir(),"results", "ABSOLUTE")),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      read.table(
        file.path(tempdir(),"results", "ABSOLUTE",
                  "DoAbsolute.YJ.ABSOLUTE.table.txt"),
        sep = "\t", header = TRUE
      ),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      list.files(file.path(tempdir(),"results", "ABSOLUTE", "summary_res")),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      read.table(
        file.path(tempdir(),"results", "ABSOLUTE", "summary_res",
                  "DoAbsolute.PP-calls_tab.txt"),
        sep = "\t", header = TRUE
      ),
      style = "serialize"
    )
    testthat::expect_snapshot_value(
      read.table(
        file.path(tempdir(),"results", "ABSOLUTE", "summary_res",
                  "DoAbsolute.PP-modes.FAILED_tab.txt"),
        sep = "\t", header = TRUE
      ),
      style = "serialize"
    )


  }
)
