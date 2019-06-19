# ROBERTS ======================================================================


#' Cluster mutations Roberts 2012
#'
#' This functions implements the method in Roberts, S. A. et al. Molecular Cell (2012) to call clusterd mutations. Regions with local higher hypermutation than expected.
#'
#' WARNING: This method considers the MR constant in each sample.
#'
#' @param vr A VRanges object
#' @param delta The maximum window size to consider mutations candidates for clustering
#' @param ce_cutoff Complex events cutoff, only last one will be kept
#' @param G_percent Percentatge of the human genome used
#' @param pval_cutoff P value to decide
#' @param dbSNP a SNPLoc object that contains the appropiate dbSNP release
#'
#' @return
#' @export
#'
#' @examples
roberts_clusters <- function(vr,
                             delta = 10000,
                             ce_cutoff = 10,
                             G_percent = 1, # 0.01 per exomes?
                             pval_cutoff = 1e-4,
                             event_categories,
                             dbSNP) {

  #browser()
  # step 0: Unique sample assumption
  stopifnot(length(unique(sampleNames(vr))) == 1)

  stopifnot(!is.unsorted(vr)) # ???

  # step 1: Total number of mutations
  n = length(vr)

  D = length(unique(sampleNames(vr))) #because we are asuming 1 sample
  G = 3234830000
  G = G * G_percent
  pi = n / (D*G)

  # step 2: filter dbSNP
  if (!is.null(dbSNP)){
    vr = roberts_filter_dbSNP(vr,dbSNP = dbSNP)
  }

  if (length(vr) == 0){
    warning("No mutations left in sample");
    return(vr)}

  # step 3: Group complex mutations
  vr = roberts_group_ce(vr = vr, ce_cutoff)

  if (length(vr) == 0){
    warning("No mutations left in sample");
    return(vr)}

  # step 4: Find recursively clusters nd compute pvalues
  # from here on important to not filter the vr! because we will lose the idx pos
  clusters = roberts_find_clust(vr = vr,delta = delta,p=pi)

  # step 5: Apply significance threshold
  idx = roberts_significance(clusts = clusters,
                                pval_cutoff = pval_cutoff )
  vr$roberts_clust = FALSE
  if (length(idx) > 0){
    vr[idx]$roberts_clust = TRUE
  }


  # call events for roberts
  pvals = as.numeric(!vr$roberts_clust) # important to reverse the mask
  events_res = detect_events(x = pvals,
                sig_cutoff = 0.5, # this val is irrelevant because !T -> 0
                event_categories = event_categories)

  vr$event_type = events_res$events
  vr$event_muts = events_res$lengths

  return(vr)
}

roberts_significance <- function(clusts,pval_cutoff){
  #browser()
  clusts = clusts[mcols(clusts)$pval < pval_cutoff]
  values_idx = purrr::map2( S4Vectors::queryHits(clusts),
                            S4Vectors::subjectHits(clusts),
                            function(from,to){
                              return(from:to)
                              })

  values_idx = unlist(values_idx)
  values_idx = values_idx[!duplicated(values_idx)]
  return(values_idx)
}

roberts_find_clust <- function(vr, delta, p) {
  clust_cand = find_pairs_VR(vr, win_length = delta)
  pvals = purrr::map2_dbl(S4Vectors::queryHits(clust_cand),
                          S4Vectors::subjectHits(clust_cand),
                          function(from, to) {
                            x = start(vr[to]) - start(vr[from]) + 1
                            # because VR 1-based
                            k = length(vr[from:to])
                            # I need the vr to be sorted!
                            pval = roberts_pvalue(x = x, k = k, p = p)
                            return(pval)
                          })

  mcols(clust_cand)$pval = pvals

  return(clust_cand)

}

roberts_group_ce <- function(vr,ce_cutoff,remove = "first"){
  pairs = find_pairs_VR(vr,ce_cutoff)

  if (length(pairs)==0){
    return(vr)
  }

  ol =length(vr)

  if (remove == "first"){
    vr = vr[-S4Vectors::queryHits(pairs)]
  } else {
    vr = vr[-S4Vectors::subjectHits(pairs)]
  }

  per = scales::percent(1-(length(vr)/ol))
  warning(glue::glue("{per} mutations removed because complex events ({ce_cutoff}bp)"))

  return(vr)
}

roberts_filter_dbSNP <- function(vr,dbSNP){
  snips = BSgenome::snpsByOverlaps(dbSNP,vr)
  ovrlaps = GenomicRanges::findOverlaps(query = snips,subject = vr)
  hitsindbSNP = S4Vectors::subjectHits(ovrlaps)
  alleles = genomicHelpersDMP::dna_codes[as.character(snips$alleles_as_ambig)]

  mask_alleles = unlist(purrr::map2(alleles,
                             hitsindbSNP,
                             function(allele,hit){
                               all(c(ref(vr[hit]),alt(vr[hit])) %in% allele)
  })) # reference and alternative have to be both in dbSNP

  hitsindbSNP = hitsindbSNP[mask_alleles]
  ol = length(vr)
  vr = vr[-hitsindbSNP,]

  per = scales::percent(1-(length(vr)/ol))
  warning(glue::glue("{per} mutations removed because present in dbSNP"))

  return(vr)
}

#' Compute the negative binomial pvalue
#'
#' @param x number of nucleotides in the cluster
#' @param k number of mutations in the cluster
#' @param p mutation probability
#'
#' @return
#' @export
#'
#' @examples
roberts_pvalue <- function(x, k, p){

  pos = x-k
  iter = 0:pos

  pvals = lapply(iter, function(j){
    choose(n = (k-1)+(j-1), k = j) * ((1-p)^j) * (p^(k-1))
  })

  res = do.call("sum",pvals)

  return(res)
}




# NIKKLIEA ======================================================================


custom_basic_clustering <- function(vr,
                                    IMD,
                                    nmuts){
  # step 0: Unique sample assumption
  stopifnot(length(unique(sampleNames(vr))) == 1)

  stopifnot(!is.unsorted(vr)) # ???

  # step 1. find pairs
  pairs = find_pairs_VR(vr,IMD)

  # step 2. get the n muts in each pair
  muts = purrr::map2_dbl(S4Vectors::queryHits(pairs),
                         S4Vectors::subjectHits(pairs),
                         function(from,to){
    muts = vr[from:to] # this needs the VR to be sorted
    return(length(muts))
  })

  # step 3. filter
  mask = muts >= nmuts
  pairs = pairs[mask]

  selected_muts = purrr::map2(S4Vectors::queryHits(pairs),
                              S4Vectors::subjectHits(pairs),
                            function(from,to){
                              return(from:to)
                            })

  selected_muts = unique(unlist(selected_muts))

  vr$custom_clust = FALSE

  if (!is.null(selected_muts)){
    vr[selected_muts]$custom_clust = TRUE
  }

  return(vr)
}


