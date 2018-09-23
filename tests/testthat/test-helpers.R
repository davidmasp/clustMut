
# helpers

context("helpers")
test_that("helpers", {
  library(VariantAnnotation)
  vr <- VRanges(seqnames = c("chr1", "chr2"),
                ranges = IRanges(c(1, 10), c(5, 20)),
                ref = c("T", "A"), alt = c("C", "T"),
                refDepth = c(5, 10), altDepth = c(7, 6),
                totalDepth = c(12, 17), sampleNames = letters[1:2],
                tumorSpecific = c(FALSE, TRUE))

  vr2 <- VRanges(seqnames = c("chr1", "chr2"),
                ranges = IRanges(c(1, 10), c(5, 20)),
                ref = c("T", "A"), alt = c("C", "T"),
                refDepth = c(5, 10), altDepth = c(7, 6),
                totalDepth = c(12, 17), sampleNames = letters[1:2],
                tumorSpecific = c(FALSE, TRUE))

  obj = list(example1 = NULL,
             example2 = vr,
             example3 = vr2)

  obj_gr = expect_warning(unlist_GR_base_list(obj))

  expect_equal(length(obj_gr),4)


})
