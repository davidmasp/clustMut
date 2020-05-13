context("test-my-test.R")

test_that("pipes works", {
  res = 2

  res %<>% magrittr::multiply_by(2)
  res2 = res %>% magrittr::multiply_by(2)
  expect_equal(res, 4)
  expect_equal(res2,8)
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

test_that("distance", {
  v1 = c(1,3,20,21,26)
  v2 = c(26,21,1,20,3)
  expect_equal(compute_m_distance(v1,k=2),c(19,18,6,18,6))
  expect_equal(compute_m_distance(v2,k=2),c(6,18,19,6,18))
})

test_that("distance tabl", {

  df = data.frame(
            V1 = c(1L, 10L, 15L, 20L, 200L),
            V2 = c(3L, 9L, 16L, 22L, 220L),
            V3 = c(2L, 8L, 18L, 17L, 215L),
            V4 = c(1L, 7L, 14L, 250L, 200L)
  )

  res = data.frame(
          V1 = c(14L, 10L, 14L, 10L, 185L),
          V2 = c(13L, 13L, 13L, 13L, 204L),
          V3 = c(15L, 10L, 10L, 15L, 198L),
          V4 = c(13L, 193L, 13L, 236L, 193L)
  )

  expect_equal(compute_distances_splited_tbl(df,factor(rep("a",5)),2),
                res)
})

test_that("distance tabl 2", {

  df = data.frame(V1 = c(10L, 20L, 10L, 20L, 25L, 1L, 21L, 210L, 200L, 15L))

  fact = c("a", "a", "b", "b", "a", "a", "b", "b", "b", "a")
  res = data.frame(
    V1 = c(10L, 10L, 11L, 180L, 10L, 14L, 11L, 189L, 180L, 10L)
  )

  expect_equal(compute_distances_splited_tbl(df,f = fact,k = 2),
               res)
})


# test_that("edit distance", {
#   str3 = c("TATATAGC",
#            "CTATATAG")
#   str2 = c(
#     "TCGTC>T",
#     "TCCTC>G",
#     "AAAAA>G"
#   )
#
#   str = c(
#     "TCGTC>T",
#     "TCCTC>C",
#     "AAAAA>G"
#   )
#
#   testthat::expect_equal(compute_edit_distance(str3)[1,2],7)
#   testthat::expect_equal(compute_edit_distance(str2)[1,2],2)
#   testthat::expect_equal(compute_edit_distance(str)[1,2],1.5)
# })


test_that("event_calling", {

  fdr_test =c(runif(10,min = 0.21),
              c(.19,.1,.05),
              runif(10,min=0.21),
              c(.1,.1,.1,.1,.1))

  categories_test = c("kataegis" = 5,"omikli" = 2)

  expected_res = list(
    c(
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      "omikli",
      "omikli",
      "omikli",
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      "kataegis",
      "kataegis",
      "kataegis",
      "kataegis",
      "kataegis"),
    c(
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      3,
      3,
      3,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      5,
      5,
      5,
      5,
      5))

  test_result =  detect_events(x = fdr_test,
                               sig_cutoff = 0.2,
                               event_categories =categories_test )

  testthat::expect_equal(expected_res[[1]],test_result[[1]] )
  testthat::expect_equal(expected_res[[2]],test_result[[2]] )

})
