

find_pairs_VR <- function(vr,win_length){
  dat_vr = vr
  win = win_length
  dat_ext = genomicHelpersDMP::extend(dat_vr,upstream = win/2,downstream = win/2)
  ovr = findOverlaps(dat_ext,dat_vr)
  pairs = ovr[!duplicated(t(apply(as.matrix(ovr), 1, sort))),]
  pairs = pairs[queryHits(pairs) != subjectHits(pairs)]
  return(pairs)
}
