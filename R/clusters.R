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
  #browser()
  library(cowplot)
  library(viridis)
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
    scale_color_viridis(option = "A",end = .9)



  return(fp)

}

compute_fdr_basic <- function(pos_distance,random_matrix){
  # this is to avoid log(0) should be done internally I guess.
  #browser()
  stopifnot(nrow(random_matrix) == length(pos_distance))

  pos_distance = pos_distance + 1
  random_matrix = random_matrix + 1

  fdr_matrix = apply(random_matrix,2,function(y){
    localFDR::compute_densityfdr(obs = log(pos_distance), null = log(y),alternative="left")
  })

  fdr_vec = apply(fdr_matrix,1,median)
  return(fdr_vec)
}



clust_dist <- function(vr,rand_df,no_cores = NULL,ce_cutoff=1){
  #browser()

  stopifnot(requireNamespace("VariantAnnotation",quietly = TRUE))

  # split the data
  split_factor = VariantAnnotation::sampleNames(vr)

  # this should generate a list of indices
  idx_split = base::split(seq_along(split_factor), split_factor)

  # this should generate a list with 2 elements divided by sample
  input_list = purrr::map(.x = idx_split,function(x){
    list(vr = vr[x],RAND = rand_df[x,])
  })

  # prepare the cluster
  if (!is.null(no_cores)){
    stopifnot(requireNamespace("parallel",quietly = TRUE))
    library(parallel)
    cl = makeCluster(no_cores)
    clusterEvalQ(cl = cl, library(clustMut))
    clusterEvalQ(cl = cl, library(VariantAnnotation))

    res = parLapply(cl = cl,
              X = input_list,
              fun = function(x){
                clust_dist_sample(vr = x$vr,rand_df = x$RAND)
              } )

    parallel::stopCluster(cl)
  } else {
    res = lapply(input_list, function(x){
      clust_dist_sample(vr = x$vr,rand_df = x$RAND)
    })
  }

  definitive = unlist_GR_base_list(res)
  return(definitive)
}







clust_dist_sample <- function(vr,rand_df,ce_cutoff = 1){
  #browser()
  stopifnot(requireNamespace("VariantAnnotation",quietly = TRUE))

  ######### ONE SAMPLE ASSUMPTION
  stopifnot(length(unique(VariantAnnotation::sampleNames(vr))) == 1 )

  #browser()

  ## FILTERS
  n_mask = grepl("N",vr$ctx)
  warning(glue::glue("{scales::percent(sum(n_mask)/length(vr))} removed due to N mask"))

  vr = vr[!n_mask]
  rand_df = rand_df[!n_mask,]

  # filter complex events

  ce_mask = mask_complex_events(vr,cutoff = ce_cutoff)

  vr = vr[!ce_mask]
  rand_df = rand_df[!ce_mask,]

  # compute distances in the random data
  split_factor = seqnames(vr)
  rand_dist = compute_distances_splited_tbl(rand_df,
                                            f = split_factor,
                                            no_cores = NULL) # explore this

  gr_dist = GenomicRanges::distanceToNearest(vr)
  # so when a chromosome only have one mutation
  vr = vr[queryHits(gr_dist)]
  rand_dist = rand_dist[queryHits(gr_dist),]
  mdist = mcols(gr_dist)$distance + 1
  random_matrix = as.matrix(rand_dist)

  fdr = compute_fdr_basic(pos_distance = mdist,
                          random_matrix = random_matrix)

  vr$fdr = fdr
  vr$dist = mdist
  vr$exp_dist = random_matrix[,sample(ncol(random_matrix),size = 1)]
  vr$tp = sum(1- vr$fdr)

  return(vr)
}

sample_free_clusters <- function(dat_gr,rand_df,plot=FALSE) {

  ######### ONE SAMPLE ASSUMPTION ###################

  stopifnot(length(unique(dat_gr$sample)) == 1 )

  #browser()
  n_mask = grepl("N",dat_gr$ctx)
  print(glue::glue("{scales::percent(sum(n_mask)/length(dat_gr))} removed due to N mask"))

  dat_gr = dat_gr[!n_mask]
  rand_df = rand_df[!n_mask,]

  # filter complex events

  ce_mask = mask_complex_events(dat_gr)

  dat_gr = dat_gr[!ce_mask]
  rand_df = rand_df[!ce_mask,]

  split_factor = seqnames(dat_gr)
  rand_dist = compute_distances_splited_tbl(rand_df,f = split_factor,no_cores = 5)

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
