# ROBERTS ======================================================================



roberts_clusters_master <- function(mode,vr,
                                    delta = 10000,
                                    ce_cutoff = 10,
                                    G_percent = 1, # 0.01 per exomes?
                                    pval_cutoff = 1e-4,
                                    event_categories,
                                    dbSNP){
  if (mode == "pairs"){
    res = roberts_clusters_pairs(vr = vr,
                                 ce_cutoff = ce_cutoff,
                                 G_percent = G_percent, # 0.01 per exomes?
                                 event_categories = event_categories,
                                 dbSNP = dbSNP)
  } else if (mode == "basic"){
    res = roberts_clusters(vr = vr,
                           delta = delta,
                           ce_cutoff = ce_cutoff,
                           G_percent = G_percent, # 0.01 per exomes?
                           pval_cutoff = pval_cutoff,
                           event_categories = event_categories,
                           dbSNP = dbSNP)
  } else {
    stop("mode not found, check again")
  }
  return(res)
}


#' Cluster binomial for pairs
#'
#'
#' This function is an adaptations of the method in Roberts, S. A. et al. Molecular Cell (2012) to call clusterd mutations based only in pairs. Regions with local higher hypermutation than expected.
#'
#' WARNING: This method considers the MR constant in each sample.
#'
#' @param vr A VRanges object
#' @param ce_cutoff Complex events cutoff, only last one will be kept
#' @param G_percent Percentatge of the human genome used
#' @param dbSNP a SNPLoc object that contains the appropiate dbSNP release
#'
#' @return
#' @export
#'
#' @examples
roberts_clusters_pairs <- function(vr,
                             ce_cutoff = 10,
                             G_percent = 1, # 0.01 per exomes?
                             event_categories,
                             dbSNP) {

  # step 0: Unique sample assumption
  stopifnot(length(unique(sampleNames(vr))) == 1)


  ol = length(vr)
  vr = removeBiallelic(vr)
  el = length(vr)
  if (ol != el){
    warning(glue::glue("Biallelic positions removed, total {scales::percent(el/ol)}"))
  }

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

  # step 4: get nearest pair
  nearest_pair = GenomicRanges::distanceToNearest(vr)
  dist = mcols(nearest_pair)$distance
  test_binom = genomicHelpersDMP::binom_test(x = rep(2,length(nearest_pair)),
                                             n = dist,
                                             p = rep(pi,length(nearest_pair)))

  test_binom$p.value.bonferroni = p.adjust(test_binom$p.value,
                                           method = "bonferroni")

  test_binom$p.value.fdr = p.adjust(test_binom$p.value,
                                    method = "fdr")

  # we set the pvalue NULL to 1, the mutations that don't have a pair
  # will have a p value == 1
  mDF = mcols(vr)
  mDF$roberts_pairs_pval = 1
  mDF$roberts_pairs_pval_bonf = 1
  mDF$roberts_pairs_pval_fdr = 1

  mDF[queryHits(nearest_pair),
      c("roberts_pairs_pval",
        "roberts_pairs_pval_bonf",
        "roberts_pairs_pval_fdr")] = test_binom[,c("p.value",
                                                   "p.value.bonferroni",
                                                   "p.value.fdr")]
  mcols(vr) = mDF

  # call events for roberts
  # bonferroni is the less inflated approach I think, the fdr is not
  # corrected enough
  events_res = detect_events(x = vr$roberts_pairs_pval_bonf,
                             sig_cutoff = 0.05,
                             event_categories = event_categories)

  vr$event_type = events_res$events
  vr$event_muts = events_res$lengths
  vr$rid = events_res$rid

  return(vr)
}


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



  # step 0: Unique sample assumption
  stopifnot(length(unique(sampleNames(vr))) == 1)

  ol = length(vr)
  vr = removeBiallelic(vr)
  el = length(vr)
  if (ol != el){
    warning(glue::glue("Biallelic positions removed, total {scales::percent(el/ol)}"))
  }

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
  # here the pairs are found as window lenght extension from each mutation.
  vr = roberts_group_ce(vr = vr, ce_cutoff * 2)

  if (length(vr) == 0){
    warning("No mutations left in sample");
    return(vr)}

  # step 4: Find recursively clusters nd compute pvalues
  # from here on important to not filter the vr! because we will lose the idx pos
  clusters = roberts_find_clust(vr = vr,delta = delta,p=pi)

  # step 5: Apply significance threshold
  idx = roberts_significance(clusts = clusters,pval_cutoff = pval_cutoff )
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
  vr$rid = events_res$rid

  return(vr)
}

roberts_significance <- function(clusts,pval_cutoff){

  clusts_lol = as.list(clusts[clusts$pval < pval_cutoff,])
  values_idx = purrr::pmap(clusts_lol,
                            function(from,to,mask,pval){
                              return(from:to)
                              })

  values_idx = unlist(values_idx)
  values_idx = values_idx[!duplicated(values_idx)]
  return(values_idx)
}

roberts_find_clust <- function(vr, delta, p) {
  clust_cand = find_all_runs(vr, IMD = delta)

  stopifnot(all(clust_cand$mask))

  lol = as.list(clust_cand)

  pvals = purrr::pmap_dbl(lol,
                          function(from, to,mask) {
                            x = start(vr[to]) - start(vr[from]) + 1
                            # because VR 1-based
                            k = length(vr[from:to])
                            # I need the vr to be sorted!
                            pval = roberts_pvalue(x = x, k = k, p = p)
                            return(pval)
                          })

  clust_cand$pval = pvals

  return(clust_cand)

}

roberts_group_ce <- function(vr,ce_cutoff,remove = "first"){
  pairs = find_all_pairs(vr,ce_cutoff)

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

  sqlS = seqlevelsStyle(vr)

  if (sqlS == "UCSC"){
    seqlevelsStyle(vr) = "NCBI"
  }

  snips = BSgenome::snpsByOverlaps(dbSNP,vr)
  ovrlaps = GenomicRanges::findOverlaps(query = snips,subject = vr)
  hitsindbSNP = S4Vectors::subjectHits(ovrlaps)
  alleles = helperMut::dna_codes[as.character(snips$alleles_as_ambig)]

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

  if (sqlS == "UCSC"){
    seqlevelsStyle(vr) = "UCSC"
  }

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

#' Compute cluster based in custom thresholds
#'
#' @param vr a VRanges object
#' @param IMD the inter-mutational distance to select mutations
#' @param nmuts the number of mutations per cluster requiered
#' @param event_categories a named vactor with the minimum number of mutations
#'
#' @return
#' @export
#'
#' @examples
custom_basic_clustering <- function(vr,
                                      IMD,
                                      nmuts,
                                      event_categories,
                                      nearest = TRUE){
  # step 0: Unique sample assumption
  stopifnot(length(unique(sampleNames(vr))) == 1)
  stopifnot(!is.unsorted(vr)) # ???

  clust_mask_rle = get_custom_cluster_rle(vr = vr,
                                          IMD = IMD,
                                          nmuts = nmuts,
                                          nearest = nearest)
  clust_mask_res = as(clust_mask_rle,"vector")

  vr$custom_clust = FALSE
  vr[clust_mask_res]$custom_clust = TRUE

  # call events for basic
  pvals = as.numeric(!vr$custom_clust) # important to reverse the mask
  events_res = detect_events(x = pvals,
                             sig_cutoff = 0.5, # this val is irrelevant because !T -> 0
                             event_categories = event_categories)

  vr$event_type = events_res$events
  vr$event_muts = events_res$lengths
  vr$rid = events_res$rid
  return(vr)

}


get_custom_cluster_rle <- function(vr,nearest,IMD,nmuts) {
  if (nearest){
    ovr = distanceToNearest(vr)

    imd_res = rep(NA,length(vr))
    imd_res[queryHits(ovr)] = mcols(ovr)$distance
  } else {
    precede_id = GenomicRanges::precede(x = vr)
    # this removes last positions in the chromosome
    rm_mask = !is.na(precede_id)
    tmp_vr = vr[rm_mask]
    precede_id = precede_id[rm_mask]
    precede_vr = vr[precede_id]

    distance(tmp_vr,precede_vr) -> imd
    imd_res = rep(NA,length(vr))
    imd_res[precede_id] = imd
  }

  clust_mask_imd <- imd_res <= IMD
  clust_mask_rle <- Rle(clust_mask_imd)
  events_mask = runLength(clust_mask_rle) >= nmuts
  runValue(clust_mask_rle) = runValue(clust_mask_rle) & events_mask

  return(clust_mask_rle)
}





removeBiallelic <- function(vr) {
  stopifnot(length(unique(as.character(sampleNames(vr)))) == 1)
  findOverlaps(vr,vr) -> test
  same_pos = queryHits(test[queryHits(test) != subjectHits(test)])
  if (length(same_pos) > 0){
    vr = vr[-same_pos]
  }
  vr
}

