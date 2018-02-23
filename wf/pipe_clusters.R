
library(clustMut)
library(parallel)

projects = fs::dir_ls(Sys.getenv("muts_path"),type = "directory")


library(progress)
pb <- progress_bar$new(
  format = "Computing [:bar] :percent in :elapsed ETA - :eta",
  total = length(projects), clear = FALSE, width= 100)

MSM_lol = purrr::map(projects,function(x){

  print(x)
  path = fs::dir_ls(x,glob="*-randomized.tsv")
  unclust_MSM = list()
  clust_MSM = list()
  #browser()

  #browser()
  # Calculate the number of cores
  no_cores <- 5
  # Initiate cluster
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, library(clustMut))


  MSM_list = parallel::parLapply(cl = cl,X = path,fun = function(y){ # check
    print(y)

    dat = readr::read_tsv(y)
    dat = parse_randommut_out(dat)

    dat_gr = dat$GR
    rand_df = dat$RAND
    tcga_name = tools::file_path_sans_ext(basename(y))
    dat_gr_cl = sample_free_clusters(dat_gr,rand_df,plot = TRUE)
    ggsave(fs::path(fs::path_dir(y),glue::glue("{tcga_name}.pdf")))

    unclust_MSM = genomicHelpersDMP::compute_MSM(dat_gr_cl[dat_gr_cl$fdr > 0.2])
    clust_MSM = genomicHelpersDMP::compute_MSM(dat_gr_cl[dat_gr_cl$fdr < 0.2])

    matrices = list(unclust = unclust_MSM,clust = clust_MSM)
    return(matrices)
  })

  stopCluster(cl)

  pb$tick()
  return(MSM_list)
})

library(purrr)

MSM_list = MSM_lol %>% unlist(recursive = F)

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


MSM_list = list(clust = as.matrix(clust_df),
                uncl = as.matrix(uncl_df))

saveRDS(MSM_list,file = "../output/msm_clusters.rds")

pca = list_matrix_pca(
  MSM_list,
  reference_matrix = "clust"
)

genomicHelpersDMP::plot_pca(pca)


