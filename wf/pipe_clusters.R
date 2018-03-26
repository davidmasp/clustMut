library(clustMut)
library(parallel)
library(magrittr)

## generate the output path
today = format(Sys.time(), "%Y%m%d")
outpath = fs::path("output",today)
fs::dir_create(outpath)
fs::dir_create(fs::path(outpath,"plots"))

path_list = fs::dir_ls(Sys.getenv("muts_path"),recursive = T,glob = "*-randomized.tsv")
#projects = projects[1:15]
print(glue::glue("Number of samples {length(path_list)}"))
unclust_MSM = list()
clust_MSM = list()

pair_selected =  c("C","G")

no_cores <- 3
# Initiate cluster
library(parallel)
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(clustMut))
clusterExport(cl = cl,varlist = "outpath")
clusterExport(cl = cl,varlist = "pair_selected")


MSM_list = parallel::parLapply(cl = cl,X = path_list,fun = function(y){ # check
#MSM_list = lapply(path_list,FUN = function(y){ # check
  #browser()
  print(y)

  dat = readr::read_tsv(y)
  dat = parse_randommut_out(dat)

  # remove the hyper-mutated
  # based in Campbell, B. B. et al. Cell doi:10.1016/j.cell.2017.09.048
  if ((length(dat$GR) / 3000) > 100) {return(NULL)}

  dat_gr = dat$GR
  rand_df = dat$RAND

  ## add filters here
  filter_mask = dat_gr$REF %in% pair_selected
  dat_gr = dat_gr[filter_mask]
  rand_df = rand_df[filter_mask,]
  ## finish filters

  tcga_name = tools::file_path_sans_ext(basename(y))
  dat_gr_cl = sample_free_clusters(dat_gr,rand_df,plot = TRUE)
  ggsave(fs::path(fs::path(outpath,"plots"),glue::glue("{tcga_name}.pdf")))
  #browser()
  if (sum(dat_gr_cl$fdr < 0.2) >= 10){
    #browser()
    psel = paste(pair_selected,collapse="")
    unclust_MSM = genomicHelpersDMP::compute_MSM(dat_gr_cl[dat_gr_cl$fdr >= 0.2 ])
    mask = colnames(unclust_MSM) %>% grepl(pattern = sprintf("[actgACTG]{1}[%s][actgACTG]{1}",psel))
    unclust_MSM = unclust_MSM[,mask,drop = FALSE]

    unclust_MSM_penta = genomicHelpersDMP::compute_MSM(dat_gr_cl[dat_gr_cl$fdr >= 0.2 ],k=2)
    mask = colnames(unclust_MSM_penta) %>% grepl(pattern = sprintf("[actgACTG]{2}[%s][actgACTG]{2}",psel))
    unclust_MSM_penta = unclust_MSM_penta[,mask,drop = FALSE]

    # the method is bi-tail so it also "detects" anti-clusters. Normally this is not an issue but just to be sure.
    # I keep it in the qqplots
    clust_MSM = genomicHelpersDMP::compute_MSM(dat_gr_cl[(dat_gr_cl$fdr < 0.2) & (dat_gr_cl$dist < mean(dat_gr_cl$dist))])
    mask = colnames(clust_MSM) %>% grepl(pattern =  sprintf("[actgACTG]{1}[%s][actgACTG]{1}",psel))
    clust_MSM =  clust_MSM[,mask,drop = FALSE]

    clust_MSM_penta = genomicHelpersDMP::compute_MSM(dat_gr_cl[(dat_gr_cl$fdr < 0.2) & (dat_gr_cl$dist < mean(dat_gr_cl$dist))],k=2)
    mask = colnames(clust_MSM_penta) %>% grepl(pattern =  sprintf("[actgACTG]{2}[%s][actgACTG]{2}",psel))
    clust_MSM_penta = clust_MSM_penta[,mask,drop = FALSE]


    matrices = list(unclust = unclust_MSM,
                    clust = clust_MSM,
                    unclust_penta = unclust_MSM_penta,
                    clust_penta = clust_MSM_penta)

    return(matrices)
  } else {
    warning(glue::glue("Sample {tcga_name} have less than 10 clustered mutations."))
    return(NULL)
  }
}) %>% purrr::keep(~!is.null(.))

stopCluster(cl)


print(glue::glue("Number of samples after cluster call {length(MSM_list)}"))
#MSM_list = MSM_lol %>% unlist(recursive = F)

MSM_list_backup = MSM_list

library(purrr)

unclust_m_list = MSM_list %>%
  purrr::map(extract2,"unclust") %>%
  purrr::map(as.data.frame)
uncl_df = do.call(rbind, unname(unclust_m_list))
colnames(uncl_df) <- glue::glue("{colnames(uncl_df)}_uncl")


clust_m_list = MSM_list %>%
  purrr::map(extract2,"clust") %>%
  purrr::map(as.data.frame)
clust_df = do.call(rbind, unname(clust_m_list))
colnames(clust_df) <- glue::glue("{colnames(clust_df)}_clust")


MSM_list_tri = list(clust = as.matrix(clust_df),
                uncl = as.matrix(uncl_df))

saveRDS(MSM_list_tri,file = fs::path(outpath,"msm_clusters_k1.rds"))

# gather penta nucl
unclust_m_list = MSM_list %>%
  purrr::map(extract2,"unclust_penta") %>%
  purrr::map(as.data.frame)
uncl_df = do.call(rbind, unname(unclust_m_list))
colnames(uncl_df) <- glue::glue("{colnames(uncl_df)}_uncl")

clust_m_list = MSM_list %>%
  purrr::map(extract2,"clust_penta") %>%
  purrr::map(as.data.frame)
clust_df = do.call(rbind, unname(clust_m_list))
colnames(clust_df) <- glue::glue("{colnames(clust_df)}_clust")


MSM_list_penta = list(clust_penta = as.matrix(clust_df),
                uncl_penta = as.matrix(uncl_df))

saveRDS(MSM_list_penta,file = fs::path(outpath,"msm_clusters_k2.rds"))

