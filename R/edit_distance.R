# EDIT distance

# str1 = c(
#   "TCGAC>T",
#   "TCCAC>G",
#   "TTTTT>G"
# )
#
# str2 = c(
#   "TCGTC>T",
#   "TCCTC>G",
#   "AAAAA>G"
# )


edit_distance_fdr <- function(vr,pairs_size, k, genome){
  #browser()
  library(VariantAnnotation)
  library(magrittr)

  K = 2*k + 1
  # one sample assumption
  stopifnot(length(unique(sampleNames(vr))) == 1)

  vr$mutid = glue::glue("{seqnames(vr)}_{start(vr)}_{alt(vr)}_{end(vr)}")


  # match the genome. This is a bit tricky, maybe should control more. see
  seqlevelsStyle(vr) = provider(genome)

  # solve this
  vr = vr[seqnames(vr) %in% genomicHelpersDMP::chromosomes_UCSC_in]

  dat_split <- base::split(vr,seqnames(vr))

  dat_split %>% purrr::map_df(function(x){ #could be paral
    #browser()
    MS = genomicHelpersDMP::get_MS_VR(x, k = k,genome = genome)
    MS = stringr::str_sub(MS,1,K)

    precomp_adist = adist(MS)
    diag(precomp_adist)=NA # this are the withitself values.

    dist= GenomicRanges::distanceToNearest(x)
    pairs = x %>% find_pairs_VR(pairs_size)

    if (length(pairs)==0){
      return(NULL)
    }

    p1_idx = queryHits(pairs)
    p2_idx = subjectHits(pairs)


    pair_edit_dist = vals = purrr::map2_dbl(
      .x = p1_idx,
      .y = p2_idx,
      function(i,j){
        precomp_adist[i,j]
    })

    #browser()
    expct = as.numeric(precomp_adist)
    expct = expct[!is.na(expct)]
    expct = sample(expct,replace = FALSE,
                   size = length(pair_edit_dist))
    # problematic line here

    mxEditDistance = max(c(expct,pair_edit_dist))

    exp = expct %>% cut(c(seq(0,
                              mxEditDistance,
                              by = 1)),
                        include.lowest = T) %>%
      table

    obs_vec = pair_edit_dist %>% cut(c(seq(0,
                              mxEditDistance,
                              by = 1)),
                        include.lowest = T)
    obs = obs_vec %>% table
    fdr_table = (exp/obs)


    fdr_raw = fdr_table[obs_vec]

    ord_st = order(pair_edit_dist)
    fdr_monotonic = pmin(cummax(fdr_raw[ord_st]),1)
    fdr_monotonic = fdr_monotonic[order(ord_st)]
    stopifnot(names(fdr_monotonic) == names(fdr_raw))

    pairs_df = pairs %>%  as.data.frame()
    pairs_df$fdr = fdr_monotonic


    pairs_df$mut1 = x[pairs_df$queryHits]$mutid
    pairs_df$mut2 = x[pairs_df$subjectHits]$mutid

    # substitute in the original version

    return(pairs_df)
  })  -> edit_fdr

  #browser()
  fdrsdf = data.frame(mutid = c(edit_fdr$mut1,edit_fdr$mut2),fdr = c(edit_fdr$fdr,edit_fdr$fdr))
  fdrsdf = fdrsdf %>%
    dplyr::group_by(mutid) %>%
    dplyr::summarise(n = n(), fdr = median(fdr))

  mcols(vr) = mcols(vr) %>% as.data.frame() %>%
    dplyr::full_join(fdrsdf,by = "mutid")

  # this is per sample so I can do this
  vr$tp = sum(1-fdrsdf$fdr,na.rm = TRUE)

  return(vr)
}


