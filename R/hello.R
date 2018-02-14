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



# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
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
  mcols(dat_gr) = meta_df

  if (nrow(rand_df) != length(dat_gr)) stop("kjdhfkhds")

  res = list(GR=dat_gr,
             RAND=rand_df)

  return(res)
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
