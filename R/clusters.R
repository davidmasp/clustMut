###############################################################################
#                         _           _
#                        | |         | |
#                     ___| |_   _ ___| |_ ___ _ __ ___
#                    / __| | | | / __| __/ _ \ '__/ __|
#                   | (__| | |_| \__ \ ||  __/ |  \__ \
#                    \___|_|\__,_|___/\__\___|_|  |___/
#
###############################################################################

plot_exp <- function(vr,filename = "expected_distances.pdf"){
  pdf(file = filename,paper = "a4r")

  vr %>% base::split(VariantAnnotation::sampleNames(.)) %>%
    purrr::walk(function(x){

      fp = plot_eedistances(x$dist,x$exp_dist,fdr = x$fdr,
                       sample = unique(VariantAnnotation::sampleNames(x)))

      print(fp)
    })

  dev.off()
}

plot_exp_list <- function(vr){

  res = vr %>% base::split(VariantAnnotation::sampleNames(.)) %>%
    purrr::map(function(x){

      fp = plot_eedistances(x$dist,x$exp_dist,fdr = x$fdr,
                            sample = unique(VariantAnnotation::sampleNames(x)))

      return(fp)
    })

  return(res)
}

plot_eedistances <- function(mdist,rdist,fdr,sample) {
  library(cowplot)
  require(ggplot2)

  nclust = sum(fdr<0.2)

  predicted_tp = round(sum(1 - fdr),digits = 2)

  clust_per = scales::percent(nclust/length(fdr))

  caption_text = glue::glue("Clust Muts = {nclust} ({clust_per})
                             {format(length(fdr)/2000,digits=2)} mutations / Mbp
                            Predicted total cluster TP {predicted_tp}")

  df = data.frame(
    fdr = fdr[order(mdist)],
    y = log10(mdist[order(mdist)]),
    x = log10(rdist[order(rdist)]))

  fp = ggplot(df,aes(x = -x,y=-y,color = fdr)) +
    geom_point() +
    geom_abline(slope = 1,linetype="dashed") +
    geom_vline(xintercept = -max(df[df$fdr < 0.2,]$x),color = "gray") +
    geom_hline(yintercept = -max(df[df$fdr < 0.2,]$y), color = "gray") +
    geom_rug(alpha = 1) +
    labs(x = "-log10(Expected distance to nearest)",
         y = "-log10(Observed distance to nearest)",
         title = sample,
         caption = caption_text)+
    viridis::scale_color_viridis(option = "A",end = .9)



  return(fp)

}

compute_fdr_basic <- function(pos_distance,random_matrix){
  # this is to avoid log(0) should be done internally I guess.
  stopifnot(nrow(random_matrix) == length(pos_distance))

  pos_distance = pos_distance + 1
  random_matrix = random_matrix + 1

  if (is.list(random_matrix)){
    stop("wrong random matrix inputed, why? line 88 clusters.R")
  }

  fdr_matrix = apply(random_matrix,2,function(y){
    compute_densityfdr(obs = log(pos_distance),
                       null = log(y),
                       alternative="left",
                       monotonic = FALSE)
  })
  fdr_vec = apply(fdr_matrix,1,median)
  return(fdr_vec)
}



#' Distance mutation clusters
#'
#' Detect clusters of mutaion based on the chromosomic distance or IMD.
#'
#' @param vr A VRanges object with multiple samples (or just one)
#' @param rand_df A randommut dataframe with the randomized values
#' @param ce_cutoff number of nucleotides that are considered to classify complex events
#' @param method fdr for local (per mutation) and FDR for tail-based (sample - wise) FDR
#' @param n Number of mutations enclosed in the IMD calculation. N = 1 for pairs of mutations.
#' @param dist_cutoff (FDR only) distance cutoff
#' @param split_factor (optional) By default it will split the VR by sample, however, to get a stratified boosting different factors can be provided. A vector equal size as the VR with the groups to split.
#'
#' @return A VR object with an extended metadata column with FDR or fdr values.
#' @export
#'
#' @examples
clust_dist <- function(vr,
                       rand_df,
                       ce_cutoff=1,
                       method="fdr", # FDR
                       n = 1,
                       dist_cutoff = NULL,
                       split_factor = NULL,
                       event_categories,
                       event_fdr){

  if (!is.null(dist_cutoff) & method == "fdr"){
    warning("A distance cutoff exist but method set to fdr, dropping distance cutoff")
    dist_cutoff = NULL
  }

  stopifnot(requireNamespace("VariantAnnotation",quietly = TRUE))

  # sort the VR
  vr_order = GenomicRanges::order(vr)
  vr = vr[vr_order]
  rand_df = rand_df[vr_order,]
  split_factor = split_factor[vr_order]

  # split the data
  if (is.null(split_factor)){
    split_factor = VariantAnnotation::sampleNames(vr)
  } else {
    split_factor = split_factor
  }


  # this should generate a list of indices
  idx_split = base::split(seq_along(split_factor), split_factor)

  # this should generate a list with 2 elements divided by sample
  input_list = purrr::map(.x = idx_split,function(x){
    list(vr = vr[x],RAND = rand_df[x,])
  })
  force(n)
  #browser()
  res = lapply(input_list, function(x){
    switch (method,
            fdr = clust_dist_sample(vr = x$vr,
                                    rand_df = x$RAND,
                                    n = n,
                                    event_categories = event_categories,
                                    event_fdr = event_fdr),
            FDR = clust_dist_sample_FDR(vr = x$vr,
                                        rand_df = x$RAND,
                                        dist_cutoff = dist_cutoff,
                                        n = n)
    )
  })

  definitive = unlist_GR_base_list(res)
  return(definitive)
}


clust_dist_sample <- function(vr,
                              rand_df,
                              ce_cutoff = 1,
                              n = 1,
                              event_categories = c("omikli" = 2,"kataegis" = 5),
                              event_fdr = 0.2){

  stopifnot(requireNamespace("VariantAnnotation",quietly = TRUE))

  ######### ONE SAMPLE ASSUMPTION
  sample_name = unique(VariantAnnotation::sampleNames(vr))
  stopifnot(length(sample_name) == 1 )

  ## FILTERS
  ## 1) N mask
  n_mask = grepl("N",vr$ctx)
  warning(glue::glue("{scales::percent(sum(n_mask)/length(vr))} removed due to N mask"))

  vr = vr[!n_mask]
  rand_df = rand_df[!n_mask,]

  if (0 == length(vr)){
    warning(glue::glue("Sample {sample_name} removed because no valid mutations were available"))
    return(NULL)
  }

  # 2) filter complex events
  ce_mask = mask_complex_events(vr, cutoff = ce_cutoff)

  vr = vr[!ce_mask]
  rand_df = rand_df[!ce_mask,]

  if (0 == length(vr)){
    warning(glue::glue("Sample {sample_name} removed because no valid mutations were available"))
    return(NULL)
  }

  # compute distances in the random data
  split_factor_sample = as.character(seqnames(vr))
  rand_dist = compute_distances_splited_tbl(rand_df,
                                            f = split_factor_sample,
                                            k = n) # explore this
  # numeric problem issue #21
  if (is.numeric(n) & n%%1==0){
    n = as.integer(n)
  }

  if (n == 1 & is.integer(n)){
    gr_dist = GenomicRanges::distanceToNearest(vr)
    # so when a chromosome only have one mutation, this should not be a problem
    # because we are removing this cases. Still good to be safe
    vr = vr[queryHits(gr_dist)]
    rand_dist = rand_dist[queryHits(gr_dist),]
    mdist = mcols(gr_dist)$distance + 1 # this comes from GR distance
    random_matrix = as.matrix(rand_dist)
  } else if (n > 1 & is.integer(n)) {
    vr = compute_distance_vr(vr = vr,enclosing = n)
    mdist = mcols(vr)$distance # here we don't need +1 because it's not GR
    random_matrix = as.matrix(rand_dist)
  } else {
    stop("n should be a positive integer number")
  }
  # so when a chromosome only have one mutation we place a NA
  # see issue #23 for more detail
  # what we do is to keep the NAs and then removing them at the fdr!

  if (!sum(is.na(mdist)) == sum(is.na(random_matrix))/ncol(random_matrix)){
    stop("Error kjsdkjhfdksa")
  }

  na_mask = !is.na(mdist)
  if (sum(na_mask) == 0){
    warning(glue::glue("Sample {sample_name} was removed because not enough mutations were found."))
    return(NULL)
  }
  mdist = mdist[na_mask]
  vr <- vr[na_mask]
  random_matrix = apply(random_matrix, 2, function(x){x[!is.na(x)]})

  if(is.null(dim(random_matrix))){
    warning(glue::glue("Sample {unique(sampleNames(vr))} discarded due to no vallid mutations."))
    return(NULL)
  }

  fdr = compute_fdr_basic(pos_distance = mdist,
                          random_matrix = random_matrix)

  MDF = mcols(vr)

  MDF$fdr = fdr
  MDF$dist = mdist
  rc= sample(ncol(random_matrix),size = 1)
  xp_dist = random_matrix[,rc]
  xp_dist = xp_dist[!is.na(xp_dist)]
  MDF$exp_dist = xp_dist
  MDF$tp = sum(1- MDF$fdr)

  events_out = detect_events(x = fdr,
                sig_cutoff = event_fdr,
                event_categories = event_categories)

  MDF$event_type = events_out$events
  MDF$event_muts = events_out$lengths
  MDF$event_rid = events_out$rid

  mcols(vr) = MDF

  return(vr)
}

sample_free_clusters <- function(dat_gr,rand_df,plot=FALSE) {

  ######### ONE SAMPLE ASSUMPTION ###################

  stopifnot(length(unique(dat_gr$sample)) == 1 )

  n_mask = grepl("N",dat_gr$ctx)
  print(
    glue::glue(
      "{scales::percent(sum(n_mask)/length(dat_gr))} removed due to N mask"))

  dat_gr = dat_gr[!n_mask]
  rand_df = rand_df[!n_mask,]

  # filter complex events

  ce_mask = mask_complex_events(dat_gr)

  dat_gr = dat_gr[!ce_mask]
  rand_df = rand_df[!ce_mask,]

  split_factor = seqnames(dat_gr)
  rand_dist = compute_distances_splited_tbl(rand_df,
                                            f = split_factor)

  gr_dist = distanceToNearest(dat_gr)
  # so when a chromosome only have one mutation
  dat_gr = dat_gr[queryHits(gr_dist)]
  rand_dist = rand_dist[queryHits(gr_dist),]
  mdist = mcols(gr_dist)$distance + 1
  random_matrix = as.matrix(rand_dist)

  fdr = compute_fdr_basic(pos_distance = mdist,
                          random_matrix = random_matrix)

  dat_gr$fdr = fdr
  dat_gr$dist = mdist

  if (plot){
    fp = plot_eedistances(mdist = dat_gr$dist,fdr = dat_gr$fdr,
                     rdist = rand_dist$R1,sample = unique(dat_gr$sample))
    print(fp)
  }

  return(dat_gr)
}


generate_mut_format <- function(gr){
  mut_code = glue::glue("{seqnames(gr)}:{start(gr)}-{end(gr)}_{gr$REF}-{gr$ALT}_{gr$sample}")
  return(mut_code)
}


write_clust_muts <- function(dat_gr,clust_mask,filename){
  gr = dat_gr[clust_mask]
  mut_code = generate_mut_format(gr)
  readr::write_lines(mut_code,path = filename)
}


#### FDR ######

clust_dist_sample_FDR <- function(vr,rand_df,ce_cutoff = 1,dist_cutoff,n=1){

  stopifnot(requireNamespace("VariantAnnotation",quietly = TRUE))

  ######### ONE SAMPLE ASSUMPTION
  stopifnot(length(unique(VariantAnnotation::sampleNames(vr))) == 1 )


  ## FILTERS
  n_mask = grepl("N",vr$ctx)
  warning(glue::glue("{scales::percent(sum(n_mask)/length(vr))} removed due to N mask"))

  vr = vr[!n_mask]
  rand_df = rand_df[!n_mask,]

  # filter complex events

  ce_mask = mask_complex_events(vr,cutoff = ce_cutoff)

  vr = vr[!ce_mask]
  rand_df = rand_df[!ce_mask,]

  if (0 == length(vr)){
    warning(glue::glue("Sample {sample_name} removed because no valid mutations were available"))
    return(NULL)
  }

  # compute distances in the random data
  split_factor = seqnames(vr)
  rand_dist = compute_distances_splited_tbl(rand_df,
                                            f = split_factor,
                                            k = n)

  # numeric problem issue #21
  if (is.numeric(n) & n%%1==0){
    n = as.integer(n)
  } else {
    stop("N needs to be numeric")
  }

  if (n == 1 & is.integer(n)){
    gr_dist = GenomicRanges::distanceToNearest(vr)
    # so when a chromosome only have one mutation
    vr = vr[queryHits(gr_dist)]
    rand_dist = rand_dist[queryHits(gr_dist),]
    mdist = mcols(gr_dist)$distance + 1 # this comes from GR distance
    random_matrix = as.matrix(rand_dist)
  } else if (n > 1 & is.integer(n)) {
    vr = compute_distance_vr(vr = vr,enclosing = n)
    mdist = mcols(vr)$distance # here we don't need +1 because it's not GR
    random_matrix = as.matrix(rand_dist)
  } else {
    stop("n should be a positive integer number")
  }
  # so when a chromosome only have one mutation we place a NA
  # see issue #23 for more detail
  # what we do is to keep the NAs and then removing them at the fdr!
  stopifnot(sum(is.na(mdist)) == sum(is.na(random_matrix))/ncol(random_matrix))
  mdist = mdist[!is.na(mdist)]
  if (sum(!is.na(mdist)) == 0){
    warning(glue::glue("Sample {sample_name} was removed because not enough mutations were found."))
    return(NULL)
  }
  random_matrix = apply(random_matrix, 2, function(x){x[!is.na(x)]})

  FDR = compute_FDR_basic(pos_distance = mdist,
                          dist_cutoff = dist_cutoff,
                          random_matrix = random_matrix)
  vr$FDR = FDR
  vr$TDIST = max(mdist[mdist <= dist_cutoff])
  return(vr)
}

compute_FDR_basic <- function(pos_distance,random_matrix,dist_cutoff){
  # check arguments
  if (is.null(dist_cutoff)){
    stop("Cutoff not provided")
  }

  if (is.null(dim(random_matrix)) || nrow(random_matrix) == 0){
    stop("Error in random matrix when computing FDR")
  }

  if (is.null(pos_distance) | nrow(random_matrix) != length(pos_distance)){
    stop("Error in the provided positions")
  }

  observed = sum(pos_distance < dist_cutoff)
  expected = apply(random_matrix, 2, function(x){sum(x < dist_cutoff)})

  expected = median(expected)

  FDR = expected / observed

  return(FDR)
}
