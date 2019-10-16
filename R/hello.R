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


# warning the context have to match the k of the one coming from RMut output
# so for now -> 3!
parse_randommut_vr <- function(dat,
                               context=NULL,
                               reference_set = c("C","A")){

  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyselect", quietly = TRUE)
  # suppressPackageStartupMessages(library(VariantAnnotation))
  # should be solved after issue #26

  if (!is.null(context)){
    rev_ctx = as.character(Biostrings::reverseComplement(DNAStringSet(reference_set)))
    ctx_in = c(make_set(x = context,simplify = T),
               make_set(x = context,simplify = T,simplify_set = rev_ctx))

    sep_pos = nchar(context) - 1
    sep = stringr::str_sub(context,start = sep_pos,sep_pos)
    mut_types = paste(dat$ctx,dat$alt,sep = sep)
    ol = nrow(dat)
    dat = dat[mut_types %in% ctx_in,]
    per = scales::percent(1 - (nrow(dat)/ol))
    warning(glue::glue("{per} of mutations removed due to context filter"))
  }

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
  ol = length(x)

  # be explicit when removing things.
  x = purrr::keep(x,~ !is.null(.))
  if (length(x) < ol){
    warning("Elements in vr list discarded due to NULL")
  }

  # mrge elements in the list.
  master_gr = do.call("c",x)

  return(master_gr)
}



VR_preprocessing <- function(file_paths,
                             pair_set,
                             alignability_mask){
  # read the filter
  if (!is.null(alignability_mask)){
    alignability_bed = rtracklayer::import.bed(alignability_mask)
    sq_st = seqlevelsStyle(alignability_bed)
  }


  library(progress)
  pb <- progress_bar$new(
    format = "Reading files :percent eta: :eta",
    total = length(file_paths), clear = FALSE)

  dat_list = purrr::map(file_paths,function(x){

    dat = readRDS(x)
    original = length(dat)

    # apparently this happens in ICGC
    mask = (ref(dat) == alt(dat))
    dat = dat[!mask]

    # check for duplicates
    mutid = paste(seqnames(dat),start(dat),sampleNames(dat),alt(dat),sep = ":")
    mask_duplicated = duplicated(mutid)
    dat = dat[!mask_duplicated]

    # remove bases outside the predefined set
    ref_in = genomicHelpersDMP::dna_codes[[pair_set]]
    mask=ref(dat) %in% ref_in
    dat = dat[mask]

    # I check if alignability filter is set
    if (!is.null(alignability_mask)){
      stopifnot(sq_st == seqlevelsStyle(dat))
      # remove chromosomes out of the bed file
      dat = dat[seqnames(dat) %in% seqnames(alignability_bed)]
      ovr = GenomicRanges::findOverlaps(query = dat,subject = alignability_bed)
      # we keep the mutations that are included in the bed file !!!
      dat = dat[S4Vectors::queryHits(ovr),]
    }

    after_filter = length(dat)
    per = (1 - (after_filter/original)) %>% scales::percent()

    print(glue::glue("{per} mutations discarded in sample {x}."))

    pb$tick()

    return(dat)
  })
  return(dat_list )
}


