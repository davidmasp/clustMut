

find_pairs_VR <- function(vr,win_length){
  dat_vr = vr
  win = win_length
  dat_ext = genomicHelpersDMP::extend(dat_vr,upstream = win/2,downstream = win/2)
  ovr = findOverlaps(dat_ext,dat_vr)
  pairs = ovr[!duplicated(t(apply(as.matrix(ovr), 1, sort))),]
  pairs = pairs[queryHits(pairs) != subjectHits(pairs)]
  return(pairs)
}


#' Find mutation pairs
#'
#' In a VR object, finds mutation pairs that enclose n number of mutations.
#'
#' @param vr a VRanges object
#' @param enclosing Number of mutations between pairs
#'
#' @return a pair object
#' @export
#'
#' @examples
find_mut_pairs_vr <- function(vr,enclosing){
  browser()
  stopifnot(length(unique(sampleNames(vr))) == 1)
  # single sample assumption
  ## optional, if you want a genomic order of the chromosomes

  # (??????????????)
  stopifnot(!is.unsorted(vr))

  ## split into a GRangesList
  ## where each element has all ranges for one chromosome
  vrl = split(vr, seqnames(vr))
  prev_idx = unlist(lapply(vrl,length))
  prev_idx = c(0,prev_idx[-length(prev_idx)])
  names(prev_idx) = names(vrl)
  prev_idx = cumsum(prev_idx)

  ## apply a function to the ranges of each chromosome
  res = lapply(names(vrl), function(x){
    browser()
    })
}



find_mut_pairs <- function(vr,enclosing,previous_idx,max_idx){
  #stopifnot(length(unique(sampleNames(vr))) == 1)
  stopifnot(length(unique(seqnames(vr))) == 1)
  # single sample assumption & chr
  if ( is.unsorted(start(vr))){warning("DATA IS UNSORTED!")}
  #vr = sort(vr)

  n=enclosing
  pairs_1 = integer(length(vr))
  pairs_2 = integer(length(vr))
  for (i in 1:length(vr)){
    pairs_1[i] = i
    if (i-n<=0){
      pairs_2[i] = i+n
    } else if (i+n > length(vr)){
      pairs_2[i] = i-n
    } else {
      #browser()
      fun = c("-","+")[nearest(x = vr[i],subject = vr[c(i-n,i+n)])]
      pairs_2[i] = do.call(fun,list(i,n))
    }
  }
  #browser()
  res = S4Vectors::SelfHits(from = pairs_1+previous_idx,
                            to = pairs_2+previous_idx,
                            nnode = max(max_idx))
  return(res)
}


find_mut_pairs_vr2 <- function(vr,enclosing){
  stopifnot(length(unique(sampleNames(vr))) == 1)
  # single sample assumption
  ## optional, if you want a genomic order of the chromosomes

  # (??????????????)
  stopifnot(!is.unsorted(vr))

  mcols(vr)$mdist = NA
  ## split into a GRangesList
  ## where each element has all ranges for one chromosome
  vrl = split(vr, seqnames(vr))

  ## apply a function to the ranges of each chromosome
  res = lapply(names(vrl), function(x){
    VR=vrl[[x]]
    stopifnot(length(unique(seqnames(VR))) == 1)
    dist = compute_m_distance(x = start(VR),k = enclosing,use = min)
    VR$mdist = dist
    return(VR)
  })

  names(res) <- NULL
  res = do.call("c",res)
  return(res)
}
