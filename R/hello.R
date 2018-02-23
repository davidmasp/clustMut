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
  colnames(meta_df) = gsub(pattern = "alt",replacement = "ALT",x = colnames(meta_df))
  colnames(meta_df) = gsub(pattern = "ref",replacement = "REF",x = colnames(meta_df))
  mcols(dat_gr) = meta_df

  if (nrow(rand_df) != length(dat_gr)) stop("kjdhfkhds")

  res = list(GR=dat_gr,
             RAND=rand_df)

  return(res)
}


mask_complex_events <- function(gr){
  #### unique sample assumtion
  stopifnot(length(unique(gr$sample)) == 1 )
  ndist = distanceToNearest(gr)
  idx = queryHits(ndist)
  ce_mask = !logical(length = length(gr)) # this also removes mutations that are unique in a particular chromosome
  ce_mask[idx] = mcols(ndist)$distance == 0
  warning(glue::glue("{scales::percent(sum(ce_mask) / length(gr))} mutations removed as complex events."))
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
  #browser()
  master_gr = x[[1]]
  for (i in 2:length(x)){
    master_gr = c(master_gr,x[[i]])
  }

  return(master_gr)
}




category_plot <- function(dat,category,values){
  library(cowplot)
  library(magrittr)
  library(purrr)
  library(dplyr)
  library(rlang)
  library(broom)
  requireNamespace("glue")


  fact = category

  dat[,fact] = as.factor(dat[,fact])
  val = values

  base_row = ggplot(data = dat,aes_string(x = fact)) +
    geom_bar(fill=NA,color="black") +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, NA),
      breaks = scales::pretty_breaks(n = 3)
    )

  boxplot = ggplot(data = dat,aes_string(x = fact,y=val)) +
    geom_violin(linetype="dashed") +
    geom_boxplot(fill = NA,color="darkred",width = 0.1) +
    theme(axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())


  lev_factor = levels(dat[,fact])

  val_factor = dat %>% split(.[,fact]) %>% map(pull,val)

  res_df = data.frame()

  done_comp = c()

  for (i in lev_factor){
    for (j in lev_factor){
      test = wilcox.test(val_factor[[i]],val_factor[[j]])
      df = broom::tidy(test)
      df$v1 = i
      df$v2 = j

      if (glue::glue("{j}vs{i}") %in% done_comp){
        next()
      } else {
        done_comp = c(done_comp,glue::glue("{i}vs{j}"))
      }

      res_df = rbind(res_df,df)
    }
  }

  res_df$v1 = factor(res_df$v1,levels = levels(dat[,fact]))
  res_df$v2 = factor(res_df$v2,levels = levels(dat[,fact]))

  plabs = symnum(
    res_df$p.value,
    corr = FALSE,
    na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", "n.s.")
  )
  plabs_legend = attr(plabs,"legend")
  res_df$pvallabs=  plabs %>% as.character()


  sigheatmap = res_df %>% ggplot(aes(
    x = v1,
    y = v2,
    fill = -log10(p.value),
    label = pvallabs
  )) +
    geom_tile(alpha=0.2,color = "black") +
    scale_fill_gradient2(high = "darkred",mid = "white",low = "blue") +
    geom_text() +
    theme(
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  theme_set(theme_cowplot(font_size = 8))

  fp = plot_grid(
    sigheatmap + labs(caption = glue::glue("Two-sided Mann-Whitney test:  {plabs_legend}")),
    boxplot,
    base_row,
    rel_heights = c(.2, .6, .2),
    ncol = 1,
    align = "v"
  )

  return(fp)

}

multi_values_plot <- function(dat,values,category){
  plot_list = list()

  for (i in values){
    #browser()
    plot_list[[i]] = category_plot(dat,
                                   category = category,
                                   values = i)
  }

  fp =  plot_grid(nrow = 1,plotlist = plot_list,labels = "AUTO")


  return(fp)


}

