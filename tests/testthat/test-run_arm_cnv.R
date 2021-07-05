testthat::test_that(
  "test run_arm_cnv()",{

    seg_cnv <- readRDS(system.file("extdata", "run_arm_cnv_example_seg_cnv.rds",
                                   package = "yjtools"))
    testthat::expect_snapshot_value(
      run_arm_cnv(seg_cnv, "CNV", "barcode"),
      style = "serialize"
    )

    testthat::expect_snapshot_value(
      run_arm_cnv(seg_cnv, "CNV"),
      style = "serialize"
    )

  }

)
