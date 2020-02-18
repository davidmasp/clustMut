

#' Find pairs within a IMD range
#'
#' @param vr a VRanges object
#' @param win_length the IMD in which pairs will be found
#'
#' @return
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
  pairs = ovr[!duplicated(t(apply(as.matrix(ovr), 1, sort))),]
  pairs = pairs[S4Vectors::queryHits(pairs) != S4Vectors::subjectHits(pairs)]
  return(pairs)
}

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
