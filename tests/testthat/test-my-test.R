context("test-my-test.R")

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("multiplication works", {
  library(VariantAnnotation)
  test_VR = VRanges(seqnames = c("chr1","chr1","chr1","chr1",
                                 "chr2","chr2","chr2","chr2"),
                    ranges =IRanges(start = c(10,12,13,20,12,14,21,22),
                                    end = c(10,12,13,20,12,14,21,22)),
                    sampleNames = c("S1","S2","S2","S1","S2","S2","S2","S3"),
                    ref = rep("A",8),
                    alt = rep("T",8)
  )
  testthat::expect_warning(mcols(distance_per_sample_VR(test_VR))$distance)
  testthat::expect_equal(mcols(distance_per_sample_VR(test_VR))$distance,
                         c(9,0,0,9,1,1,6,-1))
  })
