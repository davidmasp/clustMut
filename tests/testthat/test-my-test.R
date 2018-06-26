context("test-my-test.R")

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


# test_that("multiplication works", {
#   library(VariantAnnotation)
#   test_VR = VRanges(seqnames = c("chr1","chr1","chr1","chr1",
#                                  "chr2","chr2","chr2","chr2"),
#                     ranges =IRanges(start = c(10,12,13,20,12,14,21,22),
#                                     end = c(10,12,13,20,12,14,21,22)),
#                     sampleNames = c("S1","S2","S2","S1","S2","S2","S2","S3"),
#                     ref = rep("A",8),
#                     alt = rep("T",8)
#   )
#   testthat::expect_warning(mcols(distance_per_sample_VR(test_VR))$distance)
#   testthat::expect_equal(mcols(distance_per_sample_VR(test_VR))$distance,
#                          c(9,0,0,9,1,1,6,-1))
#   })


test_that("edit distance", {
  str3 = c("TATATAGC",
           "CTATATAG")
  str2 = c(
    "TCGTC>T",
    "TCCTC>G",
    "AAAAA>G"
  )

  str = c(
    "TCGTC>T",
    "TCCTC>C",
    "AAAAA>G"
  )

  testthat::expect_equal(compute_edit_distance(str3)[1,2],7)
  testthat::expect_equal(compute_edit_distance(str2)[1,2],2)
  testthat::expect_equal(compute_edit_distance(str)[1,2],1.5)
})
