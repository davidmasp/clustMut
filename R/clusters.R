###############################################################################
#                         _           _
#                        | |         | |
#                     ___| |_   _ ___| |_ ___ _ __ ___
#                    / __| | | | / __| __/ _ \ '__/ __|
#                   | (__| | |_| \__ \ ||  __/ |  \__ \
#                    \___|_|\__,_|___/\__\___|_|  |___/
#
###############################################################################



plot_eedistances <- function(mdist,rdist,fdr,sample) {
  #browser()
  library(cowplot)

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
         caption = glue::glue("{format(length(dat_gr)/2000,digits=2)} mutations / Mbp")) +
    scale_color_viridis_c(option = "A",end = .9)


  return(fp)

}

compute_fdr_basic <- function(pos_distance,random_matrix){
  # this is to avoid log(0) should be done internally I guess.
  pos_distance = pos_distance + 1
  random_matrix = random_matrix + 1

  fdr_matrix = apply(random_matrix,2,function(y){
    localFDR::compute_densityfdr(obs = log(pos_distance), null = log(y))
  })

  fdr_vec = apply(fdr_matrix,1,median)
  return(fdr_vec)
}


sample_free_clusters <- function(dat_gr,rand_df,plot=FALSE) {

  ######### ONE SAMPLE ASSUMPTION ###################

  stopifnot(length(unique(dat_gr$sample)) == 1 )


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




compute_fdr_parallel <- function(dat_gr,rand_dist,f,no_cores){
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("purrr", quietly = TRUE)
  requireNamespace("localFDR", quietly = TRUE)
  requireNamespace("broom", quietly = TRUE)
  #browser()
  df = data.frame(mdist = dat_gr$mdist)
  df = cbind(df,data.frame(rand_dist))
  dat_split = base::split(df,f = f)

  cl <- parallel::makeCluster(no_cores)
  parallel::clusterEvalQ(cl, {
    library(clustMut)
    library(magrittr)
    library(localFDR)
    })

  fdr_vec = parallel::parLapply(cl = cl,X = dat_split,fun = function(x){
    pos_distance = x$mdist
    random_matrix = x %>% dplyr::select(-mdist) %>% as.matrix()

    # this is to avoid log(0) should be done internally I guess.
    pos_distance = pos_distance + 1
    random_matrix = random_matrix + 1

    fdr_matrix = apply(random_matrix,2,function(y){
      localFDR::compute_densityfdr(obs = log(pos_distance), null = log(y))
    })

    fdr_vec = apply(fdr_matrix,1,median)

    return(fdr_vec)
  })

  stopCluster(cl)

  dat_gr$fdr =  unlist(fdr_vec)

  return(dat_gr)
}



compute_fdr_default <- function(dat_gr,rand_dist,f){
  requireNamespace("purrr", quietly = TRUE)
  requireNamespace("localFDR", quietly = TRUE)
  requireNamespace("broom", quietly = TRUE)
  pos_split = base::split(dat_gr$mdist,f = f)
  rand_split = split(rand_dist,f=f)
  #browser()
  lol = list(x = pos_split,y = rand_split,k = names(pos_split), z = names(rand_split))
  fdr_vec = purrr::pmap(.l = lol,.f = function(x,y,k,z){
    if(k != z){stop("lists do not match")}

    pos_distance = x
    random_matrix = as.matrix(y)

    # this is to avoid log(0) should be done internally I guess.
    pos_distance = pos_distance + 1
    random_matrix = random_matrix + 1

    fdr_matrix = apply(random_matrix,2,function(y){
      localFDR::compute_densityfdr(obs = log(pos_distance), null = log(y))
    })

    fdr_vec = apply(fdr_matrix,1,median)

    return(fdr_vec)

  } )

  dat_gr$fdr =  unlist(fdr_vec)

  return(dat_gr)
}

compute_fdr <- function(dat_gr,rand_dist,f,no_cores=NULL){

  if (is.list(f)) f <- interaction(f, drop = TRUE)

  if (!is.null(no_cores)){
    dat_fdr = compute_fdr_parallel(dat_gr,rand_dist,f,no_cores = no_cores)
  } else {
    dat_fdr = compute_fdr_default(dat_gr,rand_dist,f)
  }

  return(dat_fdr)

}



write_clust_muts <- function(dat_gr,clust_mask,filename){
  gr = dat_gr[clust_mask]

  mut_code = glue::glue("{seqnames(gr)}:{start(gr)}-{end(gr)}_{gr$REF}-{gr$ALT}_{gr$sample}")

  readr::write_lines(mut_code,path = filename)
}
