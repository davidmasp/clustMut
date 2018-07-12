###############################################################################
#                 _          _
#                | |        | |
#                | |__   ___| |_ __   ___ _ __ ___
#                | '_ \ / _ \ | '_ \ / _ \ '__/ __|
#                | | | |  __/ | |_) |  __/ |  \__ \
#                |_| |_|\___|_| .__/ \___|_|  |___/
#                              | |
#                              |_|
#
###############################################################################


parse_randommut_vr <- function(dat){
  #browser()
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyselect", quietly = TRUE)
  # suppressPackageStartupMessages(library(VariantAnnotation))
  # should be solved after issue #26

  dat_vr = VRanges(seqnames = dat$chr,
                   ranges = IRanges(start = dat$end,
                                    end = dat$end),
                   ref = dat$ref,
                   alt = dat$alt,
                   sampleNames = dat$sample
                   )

  dat_vr$ctx = dat$ctx

  rand_df = dat %>% dplyr::select(tidyselect::matches("R[0-9]+"))
  meta_df = dat %>% dplyr::select(ctx)
  mcols(dat_vr) = meta_df

  if (nrow(rand_df) != length(dat_vr)) stop("kjdhfkhds")


  res = list(VR=dat_vr,
             RAND=rand_df)

  return(res)
}



# dat = readr::read_tsv(path)
parse_randommut_out <- function(dat){
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyselect", quietly = TRUE)
  library(GenomicRanges)
  dat_gr = GRanges(seqnames = dat$chr,
                   ranges = IRanges(start = dat$end,
                                    end = dat$end),
                   strand = ifelse(grepl(pattern = "[1+*]",x = dat$strand),
                                   "+",
                                   "-"))
  rand_df = dat %>% dplyr::select(tidyselect::matches("R[0-9]+"))
  meta_df = dat %>% dplyr::select(-chr,-start,-end,-strand,-tidyselect::matches("R[0-9]+"))
  colnames(meta_df) = gsub(pattern = "alt",replacement = "ALT",x = colnames(meta_df))
  colnames(meta_df) = gsub(pattern = "ref",replacement = "REF",x = colnames(meta_df))
  mcols(dat_gr) = meta_df

  if (nrow(rand_df) != length(dat_gr)) stop("kjdhfkhds")

  res = list(GR=dat_gr,
             RAND=rand_df)

  return(res)
}


#' Mask complex events and unique mutations in a chromosome
#'
#' @param gr a VRanges object (unisample)
#' @param cutoff minimum distance to consider as a complex event
#'
#' @return
#' @export
#'
#' @examples
mask_complex_events <- function(gr, cutoff = 1) {
  #### unique sample assumtion
  if (is(gr, "VRanges")) {
    stopifnot(length(unique(VariantAnnotation::sampleNames(gr))) == 1)

  } else if ("sample" %in% colnames(mcols(gr))) {
    stopifnot(length(unique(gr$sample)) == 1)
  } else {
    warning("Format not identified, check samples")
  }

  ndist = distanceToNearest(gr)
  idx = queryHits(ndist)
  # this also removes mutations that are unique in a particular chromosome
  ce_mask = !logical(length = length(gr))
  ce_mask[idx] = mcols(ndist)$distance < cutoff
  warning(
    glue::glue(
      "{scales::percent(sum(ce_mask) / length(gr))} mutations removed as complex events."
    )
  )
  return(ce_mask)
}


binom_test <- function(x,n,p,...){
  if (!all(c(requireNamespace("broom", quietly = TRUE),
             requireNamespace("purrr", quietly = TRUE)))){
    print("Dependencies failed.")
  }
  res = purrr::map_df(1:length(x),function(i){
    broom::tidy(binom.test(x = x[i],n = n[i],p = p[i],...))
  })

  return(res)
}


unlist_GR_base_list <- function(x){
  names(x) = NULL
  master_gr = do.call("c",x)

  return(master_gr)
}



