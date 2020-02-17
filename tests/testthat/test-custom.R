#!/usr/bin/env Rscript

# Name: test custom methods
# Author: DMP
# Description: Test method for custom cluster detection

# custom -----------------------------------------------------------------


context("custom cluster detection")
test_that("custom", {
  library(VariantAnnotation)
  dat_vr = VRanges(
    seqnames = c(1,1,1,1,1,1,1,2,2),
    ranges = IRanges(
      start = c(100,150,200,250,350,500,505,100,500),
      width = 1
    ),
    ref = c("C"),
    alt = "A"
  )

  custom_basic_clustering(vr = dat_vr,
                            IMD = 50,
                            nmuts = 4,
                            event_categories = c("clust" = 3),
                            nearest = TRUE) -> tmp_res

  expect_equal(sum(tmp_res$custom_clust),4)

  custom_basic_clustering(vr = dat_vr,
                            IMD = 50,
                            nmuts = 2,
                            event_categories = c("clust" = 3),
                            nearest = TRUE) -> tmp_res

  expect_equal(sum(tmp_res$custom_clust),6)

  custom_basic_clustering(vr = dat_vr,
                            IMD = 50,
                            nmuts = 2,
                            event_categories = c("clust" = 3),
                            nearest = FALSE) -> tmp_res

  expect_equal(sum(tmp_res$custom_clust),3)


})
