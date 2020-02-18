

#' Find all pairs within a window length
#'
#' It finds all the possible pairs of mutations which are found within a
#' extended window from each mutation. The total window length
#' is extended half upstream and half downstream.
#'
#' @param vr a VRanges object
#' @param win_length the total window length where to find pairs
#'
#' @return
#'
#' Returns a overlaps object with all the possible pairs. They are
#' de-duplicated meaning that each lines represents only one pair of
#' mutations.
#'
#' @export
#'
#' @examples
find_all_pairs <- function(vr,win_length){
  dat_vr = vr
  win = win_length
  dat_ext = genomicHelpersDMP::extend(dat_vr,
                                      upstream = win/2,
                                      downstream = win/2)
  ovr = GenomicRanges::findOverlaps(dat_ext,
                                    dat_vr)
  dup_mask = as.vector(!duplicated(t(apply(as.matrix(ovr), 1, sort))))
  pairs = ovr[dup_mask]
  pairs = pairs[S4Vectors::queryHits(pairs) != S4Vectors::subjectHits(pairs)]
  return(pairs)
}

#' Find all mutation runs
#'
#' It search all the possible run-pairs according to a minimum inter-mutational
#' distance (IMD) in the vr objects and selects the ones
#' which have at least a given number of mutations.
#'
#' A run-pair is a set of consequtive pairs which each consists in less than
#' a given IMD.
#'
#' @param vr a VRanges object
#' @param IMD The IMD threshold (equal or smaller)
#'
#' @return
#'
#' A data frame with 3 columns, a from value which returns the original idx
#' of the run-pair, a to column with the final idx of the run and a all TRUE
#' mask column which is used for validation.
#'
#' @export
#'
#' @examples
find_all_runs <- function(vr,IMD){

  get_custom_cluster_rle(vr = vr,
                         IMD = IMD,
                         nmuts = 2,
                         nearest = TRUE) -> res_rle

  stopifnot(sum(runLength(res_rle)) == length(vr))

  pairs_df = data.frame(
    from = S4Vectors::start(res_rle),
    to = S4Vectors::end(res_rle),
    mask = S4Vectors::runValue(res_rle)
  )

  return(pairs_df[pairs_df$mask,])
}
