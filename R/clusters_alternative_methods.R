# roberts


roberts_clusters <- function(vr,
                             delta = 10000,
                             ce_cutoff = 10,
                             G_percent = 1, # 0.01 per exomes?
                             pval_cutoff = 1e-4,
                             dbSNP) {

    # step 0: Unique sample assumption
  stopifnot(length(unique(sampleNames(vr))) == 1)

  # step 1: Total number of mutations
  n = length(vr)

  D = length(unique(sampleNames(vr))) #because we are asuming 1 sample
  G = 3234830000
  G = G * G_percent
  pi = n / (D*G)

  # step 2: filter dbSNP
  vr = roberts_filter_dbSNP(vr,dbSNP = dbSNP)

  # step 3: Group complex mutations
  vr = roberts_group_ce(vr = vr, ce_cutoff)

  # step 4: Find recursively clusters nd compute pvalues
  # from here on important to not filter the vr! because we will lose the idx pos
  clusters = roberts_find_clust(vr = vr,delta = delta,p=pi)

  # step 5: Apply significance threshold
  idx = roberts_significance(clusts = clusters,
                                pval_cutoff = pval_cutoff )
  vr$roberts_clust = FALSE
  vr[idx]$roberts_clust = TRUE

  return(vr)
}

roberts_significance <- function(clusts,pval_cutoff){
  clusts = clusts[mcols(clusts)$pval < pval_cutoff]
  values_idx = purrr::map2( queryHits(clusts),subjectHits(clusts),
                            function(from,to){
                              return(from:to)
                              })

  values_idx = unlist(values_idx)
  values_idx = values_idx[duplicated(values_idx)]
  return(values_idx)
}

roberts_find_clust <- function(vr,delta,p){
  clust_cand = find_pairs_VR(vr,win_length = delta)
  pvals = purrr::map2_dbl(queryHits(clust_cand),
              subjectHits(clust_cand),
              function(from,to){
    x = start(vr[to]) - start(vr[from]) + 1 # because VR 1-based
    k = length( vr[from:to])
    pval = roberts_pvalue(x = x,k = k,p = p)
    return(pval)
  })

  mcols(clust_cand)$pval = pvals

  return(clust_cand)

}

roberts_group_ce <- function(vr,ce_cutoff,remove = "first"){
  pairs = find_pairs_VR(vr,ce_cutoff)

  ol =length(vr)

  if (remove == "first"){
    vr = vr[-queryHits(pairs)]
  } else {
    vr = vr[-subjectHits(pairs)]
  }

  per = scales::percent(1-(length(vr)/ol))
  warning(glue::glue("{per} mutations removed because complex events ({ce_cutoff}bp)"))

  return(vr)
}

roberts_filter_dbSNP <- function(vr,dbSNP){
  snips = BSgenome::snpsByOverlaps(dbSNP,vr)
  ovrlaps = findOverlaps(query = snips,subject = vr)
  hitsindbSNP = subjectHits(ovrlaps)
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

#' Title
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
  #browser()
  pos = x-k
  iter = 0:pos

  pvals = lapply(iter, function(j){
    choose(n = (k-1)+(j-1), k = j) * ((1-p)^j) * (p^(k-1))
  })

  res = do.call("sum",pvals)

  return(res)
}
