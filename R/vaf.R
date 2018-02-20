
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
# this functions are used in the clean.R script, don't remove

#' Compute VAF
#'
#' @param alt_reads reads of the alternative allele, note that this definition may mean different things.
#' @param total_reads total number of reads
#' @param method a character string indicating which VAF method should be used. One of 'symmetric', 'reduction' or 'simple' (default)
#'
#' @return
#' @export
#'
#' @examples
compute_VAF <- function(alt_reads,total_reads,method = "simple"){
  VAF_vec = alt_reads/total_reads
  if(method == "symmetric"){
    VAF_vec[VAF_vec > 0.5] = 1 -VAF_vec[VAF_vec>0.5]
  } else if (method == "reduction") {
    VAF_vec[VAF_vec > .5] = ">.5" # SELF NOTE! improve this #CAREFULL
  } else if (method =="simple"){
    VAF_vec = VAF_vec
  } else {
    stop("Method not yet implemented. Available methods 'symmetric', 'reduction' and 'simple'")
  }
}


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
                           unclust_distance) {

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
        diff[diff < unclust_distance] = NA

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
