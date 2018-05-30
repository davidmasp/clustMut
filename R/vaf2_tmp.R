# library(magrittr)
# library(VariantAnnotation)
#
# simulation_size_input = 5000
# dat = genomicHelpersDMP::parse_strelka_VR(dat)
#
# seqlevelsStyle(dat) <- "UCSC"
#
#
# idx = queryHits(findOverlaps(query = dat,subject = alignability))
#
# stopifnot(length(idx) == length(unique(idx))) # safety check
#
# (length(idx) / length(dat)) %>% scales::percent() %>% glue::glue("{.} included")
#
# dat = dat[idx]
#
# neighbor = distanceToNearest(dat)
# ce_mask = queryHits(neighbor[mcols(neighbor)$distance > 1])
#
#
# (1 - (length(ce_mask) / length(dat))) %>% scales::percent() %>% glue::glue("{.} removed due to C.E.")
#
# dat = dat[ce_mask]
#
# genomicHelpersDMP::plot_rainPlot(dat)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ #

local_vaflr_fdr <- function(vr, simulation_size_input) {

  dat = vr

  dat$mutid = glue::glue("{seqnames(dat)}_{start(dat)}_{alt(dat)}_{end(dat)}")
  dat$RL = dat$RR + dat$AR

  dat_split = dat %>% base::split(seqnames(dat))



  source("R/pairs.R")

  nclust = dat_split %>% purrr::map_df(function(dat){
    #browser()


    simulation_size = min(simulation_size_input,length(dat))

    vaf1 = sample(dat$VAF,size = simulation_size)
    vaf2 = sample(dat$VAF,size = simulation_size)

    dneg = abs(log2(vaf1/vaf2))


    pairs = dat %>% find_pairs_VR(10000)

    #print(seqnames(dat))
    #browser()

    if (length(pairs)==0){
      return(NULL)
    }

    vaf1 = dat[queryHits(pairs)]$VAF
    vaf2 = dat[subjectHits(pairs)]$VAF

    RL1 = dat[queryHits(pairs)]$RL
    RL2 = dat[subjectHits(pairs)]$RL
    #browser()
    dist = start(dat[subjectHits(pairs)]) - start(dat[queryHits(pairs)])

    dpos = 1:length(vaf1) %>% purrr::map(function(i){

      #browser()
      vaf_average = mean(vaf1[i],vaf2[i])
      shared_reads = max(((50 - dist[i])/(50 + dist[i])),0) # careful here!

      shared_size = round(simulation_size * shared_reads)
      same_reads = double(shared_size)
      new_simulation_size = simulation_size - shared_size
      #browser()
      vaf_simulated_1 = rbinom(n = new_simulation_size,
                               size = RL1[i],
                               prob = vaf_average) / RL1[i]
      vaf_simulated_2 = rbinom(n = new_simulation_size,
                               size = RL2[i],
                               prob = vaf_average) / RL2[i]



      vaf_simulated_LR = abs(log2(vaf_simulated_1/vaf_simulated_2))

      res = c(vaf_simulated_LR,same_reads)


      return(res)
    })

    neg = dneg %>% cut(c(seq(0,0.9,by = 0.1),Inf),include.lowest = T) %>% table
    vaflr = abs(log2(vaf1/vaf2))
    bins = vaflr %>% cut(c(seq(0,0.9,by = 0.1),Inf),include.lowest = T)

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
  fdrsdf = fdrsdf %>% dplyr::group_by(mutid) %>% dplyr::summarise(n = n(), vaf_fdr_median = median(fdr))


  mcols(dat) = mcols(dat) %>% as.data.frame() %>% dplyr::full_join(fdrsdf)
  return(dat)

}

