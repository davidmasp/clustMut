
################################################################################
#                   __      __     ______
#                   \ \    / /\   |  ____|
#                    \ \  / /  \  | |__
#                     \ \/ / /\ \ |  __|
#                      \  / ____ \| |
#                       \/_/    \_\_|
#
################################################################################


## VAF general ====================


compute_vafLR <- function(vaf1,vaf2){
  vaflr = abs(log(vaf1 / vaf2,
                  base = 2))
  return(vaflr)
}


## VAF generation (Positive) ====================


compute_shared_reads <- function(mdist,rl=50){
  same_read = (rl - mdist)/(rl + mdist)
  same_read[same_read<0] = 0
  return(same_read)
}

vaf_generation_binomial  <- function(ar, rr, mdist) {
  tr = ar + rr
  vaf = ar/tr

  ar_prime = rbinom(n = length(ar),size = tr,prob = vaf)

  p_shared = compute_shared_reads(mdist)

  vaf_prime = ((p_shared * ar)/ tr) + (((1 - p_shared)*ar_prime)/ tr)
}

## VAFlr random distribution  =============
# CAREFUL SAMPLES NOT CON SIDERED!!!!

vaflr_random <- function(vaf_vector){
  vaf1 = vaf_vector
  vaf2 = sample(vaf_vector)

  vaflr = compute_vafLR(vaf1,vaf2)

  #vaflr_df = bin_vaflr(vaflr,groups = "random")

  return(vaflr)
}


## VAFlr pipes  =============

vaflr_positive <- function(dat_gr_dist,clust_mask) {
  #browser()
  dat_gr_cl = dat_gr_dist[clust_mask]
  vaf_prime = vaf_generation_binomial(ar = dat_gr_cl$AR,
                              rr = dat_gr_cl$RR,
                              mdist = dat_gr_cl$distance)
  vaflr = compute_vafLR(vaf1 = dat_gr_cl$VAF,vaf2 =vaf_prime)
  return(vaflr)
}

vaflr_negative <- function(clust_mask,
                           f,
                           gr,
                           unclust_distance,
                           max_distance=1000000) {

  #dat_gr_splited = split(gr,f)
  cmask_splited = base::split(clust_mask,f)
  pos_splited = base::split(end(gr),f)
  vaf_splited = base::split(gr$VAF,f)

  lol = list(cmask_splited,
             pos_splited,
             vaf_splited,
             names(cmask_splited),
             names(vaf_splited),
             names(pos_splited))

  # very critical function WRITE TEST!
  pp = purrr::pmap(.l = lol,.f=function(x,y,z,q,k,g){
    # x cmask
    # y pos
    # z vaf

    if(any(k!=c(q,g))) stop("kjdkfjs")
    if (!any(x)){
      return(NULL) # no clusters found
    } else if (length(x) == sum(x)){
      warning(glue::glue("not cls muts found in {q}"))
      return(NULL)
    } else {
      idx = which(x)
      vaf_cluster = z[idx]

      vaf_unclust = double(length = length(idx))

      for (idx_out in seq_along(idx)){
        #browser()
        i = idx[idx_out]
        diff = abs(y - y[i])
        diff[diff < unclust_distance | diff > max_distance] = NA
        j = which.min(diff)

        #browser()
        #print(j)
        vaf_unclust[idx_out] = z[j]
      }
      #browser()
      return(list(cluster = vaf_cluster,
                  unclust = vaf_unclust))
    }
  }) %>% purrr::keep(~!is.null(.))


  vaf_clust = pp %>% purrr::map(extract2,"cluster") %>% unlist()
  vaf_unclust = pp %>% purrr::map(extract2,"unclust") %>% unlist()

  vaflr = compute_vafLR(vaf1 = vaf_clust,vaf2 = vaf_unclust)
  return(vaflr)

}



vaflr_observed <- function(clust_mask,
                           f,
                           gr) {
  cluster_split = base::split(clust_mask,f)
  pos_split = base::split(end(gr),f)
  mdist_split = base::split(gr$distance,f)
  vaf_split = base::split(gr$VAF,f)

  lol = list(cluster_split,
             pos_split,
             mdist_split,
             vaf_split)

  pp = purrr::pmap(.l = lol,.f = function(cl,pos,mdist,vaf){

    if(sum(cl)<1){return(NULL)}

    idx = which(cl)

    vaf_out_first = vaf[idx]

    #browser()
    vaf_out_last = idx %>% purrr::map_dbl(function(i){

      npos = pos[i] + mdist[i]
      dpos = pos[i] - mdist[i]
      pair_idx = which(pos %in% c(npos,dpos))

      # can happen
      if (length(pair_idx) == 2) {
        warning("symetric cluster, chosing one randomly")
        pair_idx = sample(x = pair_idx,size = 1)
      }
      if (length(vaf) < pair_idx){
        browser()
        stop("jhdfksfu")
        }
      return( vaf[pair_idx] )
    })
    return(data.frame(first = vaf_out_first,
                      last = vaf_out_last))
  }) %>% purrr::keep(~!is.null(.))

  vaf_obs = dplyr::bind_rows(pp)

  vaflr = compute_vafLR(vaf1 = vaf_obs$first,vaf2 = vaf_obs$last)
  return(vaflr)
}


load_vaf_fdr <- function(dat_gr,
                         clust_mask,
                         break_points) {
  #browser()

  vaflr = vaflr_positive(dat_gr_dist = dat_gr,clust_mask = clust_mask)
  vaflr = vaflr[!is.na(vaflr)] # check where these come from
  vaflr = vaflr[!is.infinite(vaflr)] # check where these come from
  #' We can see an enrichment of the low vaflr section.
  vaflr %>% density() %>% plot
  vaflr %>% hist(breaks = 30)

  #' We transform the values using the break_points previously defined
  positive_df = data.frame(vaflr = vaflr,
                           case = "positive")
  positive_df$cutpoint = cut(positive_df$vaflr,breaks = break_points,include.lowest = T)




  # negative ==========================================


  #### ==== function starts here
  ###  cluster mask
  ###  position
  ###  factor to split
  ###  vaf


  f = interaction(list(dat_gr$sample)) # here sample should also
  vaflr = vaflr_negative(clust_mask =clust_mask,
                         f = f,
                         gr = dat_gr,
                         unclust_distance = 10000 )
  #' We can see an enrichment of the low vaflr section.
  vaflr %>% abs() %>%  density() %>% plot
  vaflr %>% hist(breaks = 30)


  #' We transform the values using the break_points previously defined
  negative_df = data.frame(vaflr = vaflr,
                           case = "negative")
  negative_df$cutpoint = cut(negative_df$vaflr,breaks = break_points,include.lowest = T)


  # kernel estimator
  max_lr = max(max(negative_df$vaflr,positive_df$vaflr))

  pos_den = density(positive_df$vaflr,from=0,to=max_lr)
  neg_den = density(negative_df$vaflr,from=0,to=max_lr)

  fdr = neg_den$y /  (neg_den$y + pos_den$y)

  fdr_df = data.frame(fdr = fdr,
                      vaflr = neg_den$x)

  p3 = ggplot(fdr_df,aes(x = vaflr,y = fdr)) + geom_line() +
    geom_vline(xintercept = 0.5,linetype="dashed",color="gray") +
    ylim(0,1)

  # plots ==============

  df = rbind(positive_df,negative_df)

  p1 = df %>% ggplot(aes(x = vaflr,color = case)) + geom_density() +
    theme(legend.position = "top") +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = 0.5,linetype="dashed",color="gray")

  library(cowplot)
  library(dplyr)
  vaf_df = df %>% group_by(case) %>% mutate(total = n()) %>%  ungroup() %>%
    dplyr::group_by(cutpoint,case) %>%
    dplyr::summarise(n = n(),
                     per = n/unique(total)) %>% ungroup()
  p0 = vaf_df %>% ggplot(aes(x = cutpoint,y = n, fill = case)) +
    geom_bar(stat = "identity") +
    theme(legend.position = "top") +
    scale_fill_brewer(palette = "Set1")

  vaf_fdr = vaf_df %>% select(case,n,cutpoint) %>% tidyr::spread(key = "case",value ="n")

  vaf_fdr$Obs =  (vaf_fdr$positive + vaf_fdr$negative)
  vaf_fdr$fdr = vaf_fdr$negative / vaf_fdr$Obs

  binom_df = binom_test(x = vaf_fdr$negative,
                        n = vaf_fdr$Obs,
                        p = vaf_fdr$fdr)
  vaf_fdr = cbind(binom_df,vaf_fdr)

  cutpoints = vaf_fdr$cutpoint %>% stringr::str_split(",") %>%
    purrr::map_df(function(x){
      res = stringr::str_extract(string = x,pattern = "[:digit:][.]*[:digit:]*")
      data.frame(matrix(as.numeric(res),ncol = 2),stringsAsFactors = F)
    })

  colnames(cutpoints) <- c("lower","upper")

  vaf_fdr = cbind(cutpoints,vaf_fdr)

  p2 = vaf_fdr %>% ggplot(aes(x = cutpoint,y = fdr)) +
    geom_point() +
    geom_errorbar(aes(ymin=conf.low,ymax=conf.high),width = 0.2) +
    geom_path(group=1) +
    ylim(0,1)

  fp = plot_grid(p0,p1,p2,p3,ncol = 2,labels = "AUTO",align="v")

  bin_fdr = vaf_fdr
  den_fdr = fdr_df
  return(list(bin_fdr = bin_fdr,den_fdr = den_fdr,plot = fp))
}


### testing and bining =====================


bin_vaflr <- function(vaflr,
                      groups,
                      breaks = c(0,0.1,.2,.3,.4,.5,Inf)){

  require(purrr)

  df = data.frame(vaflr = vaflr,
                  group = groups)

  res = df %>% split(.$group) %>%
    map_df(function(x){
      res_df = cut(x$vaflr,breaks = breaks,include.lowest = TRUE) %>%
        table %>%
        as.data.frame()
      colnames(res_df) <- c("Intervals","N")
      if (length(unique(x$group))>1) stop("kjshfkjdh")
      res_df$group = unique(x$group)
      return(res_df)
    })
  return(res)
}


lm_bin_vaflr <- function(vaflr,
                         groups,
                         breaks = c(0,0.1,.2,.3,.4,.5,Inf),
                         intercept = TRUE){
  df = bin_vaflr(vaflr,
                 groups,
                 breaks = breaks)

  if (intercept){
    res = df %>% spread("group","N") %>% lm(formula = observed ~ negative + positive)
  } else {
    res = df %>% spread("group","N") %>% lm(formula = observed ~ negative + positive + 0)
  }


  return(res)
}



# VAF per mutation ======================================


#' Detect clusters using per mutation VAF fdr
#'
#' @param vr A VRanges object with somatic mutations
#' @param simulation_size_input The number of simulated VAFs to generate for the expected and idealized distribution
#' @param pairs_size The lenght (in bp) used to find possible pairs of mutations.
#'
#' @return A VRanges object with added columns (n and vaf_fdr_median)
#'
#' @export
#'
#' @examples
local_vaflr_fdr <- function(vr,
                            simulation_size_input,
                            pairs_size = 10000) {

  stopifnot(requireNamespace("VariantAnnotation",quietly = TRUE))
  stopifnot(requireNamespace("purrr",quietly = TRUE))
  stopifnot(requireNamespace("dplyr",quietly = TRUE))

  stopifnot(is(vr,"VRanges"))

  dat = vr

  dat$mutid = glue::glue("{seqnames(dat)}_{start(dat)}_{alt(dat)}_{end(dat)}")
  dat$RL = dat$RR + dat$AR

  dat_split = dat %>% base::split(seqnames(dat))

  nclust = dat_split %>% purrr::map_df(function(dat){
    simulation_size = min(simulation_size_input,length(dat))

    vaf1 = sample(dat$VAF,size = simulation_size)
    vaf2 = sample(dat$VAF,size = simulation_size)

    dneg = abs(log2(vaf1/vaf2))

    pairs = dat %>% find_pairs_VR(pairs_size)

    if (length(pairs)==0){
      return(NULL)
    }

    vaf1 = dat[queryHits(pairs)]$VAF
    vaf2 = dat[subjectHits(pairs)]$VAF

    RL1 = dat[queryHits(pairs)]$RL
    RL2 = dat[subjectHits(pairs)]$RL

    dist = start(dat[subjectHits(pairs)]) - start(dat[queryHits(pairs)])

    dpos = 1:length(vaf1) %>% purrr::map(function(i){

      vaf_average = mean(vaf1[i],vaf2[i])
      shared_reads = max(((50 - dist[i])/(50 + dist[i])),0) # careful here!

      shared_size = round(simulation_size * shared_reads)
      same_reads = double(shared_size)
      new_simulation_size = simulation_size - shared_size

      vaf_simulated_1 = stats::rbinom(n = new_simulation_size,
                                      size = RL1[i],
                                      prob = vaf_average) / RL1[i]
      vaf_simulated_2 = stats::rbinom(n = new_simulation_size,
                                      size = RL2[i],
                                      prob = vaf_average) / RL2[i]

      vaf_simulated_LR = abs(log2(vaf_simulated_1/vaf_simulated_2))

      res = c(vaf_simulated_LR,same_reads)

      return(res)
    })

    neg = dneg %>% cut(c(seq(0,0.9,by = 0.1),Inf),include.lowest = T) %>%
      table

    vaflr = abs(log2(vaf1/vaf2))

    bins = vaflr %>%
      cut(c(seq(0,0.9,by = 0.1),Inf),include.lowest = T)

    stopifnot(length(bins) == length(dpos))

    fdrs = 1:length(bins) %>% purrr::map_dbl(function(i){

      pos = dpos[[i]] %>% cut(c(seq(0,0.9,by = 0.1),Inf),include.lowest = T) %>% table

      fdr = neg[bins[i]]/(pos[bins[i]] + neg[bins[i]])
    })

    pairs_df = as.data.frame(pairs)
    pairs_df$mut1 = dat[pairs_df$queryHits]$mutid
    pairs_df$mut2 = dat[pairs_df$subjectHits]$mutid
    pairs_df$fdr = fdrs
    pairs_df$dist = dist
    return(pairs_df)
  })

  fdrsdf = data.frame(mutid = c(nclust$mut1,nclust$mut2),fdr = c(nclust$fdr,nclust$fdr))
  fdrsdf = fdrsdf %>%
    dplyr::group_by(mutid) %>%
    dplyr::summarise(n = n(), vaf_fdr_median = median(fdr))


  mcols(dat) = mcols(dat) %>% as.data.frame() %>% dplyr::full_join(fdrsdf)
  return(dat)

}

