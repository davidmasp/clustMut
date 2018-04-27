source("wf/snip/load_cmc.R")

pca = genomicHelpersDMP::list_matrix_pca(
  MSM_list,
  reference_matrix = "clust",
  scale = FALSE
)

pca_raw = genomicHelpersDMP::plot_pca(pca,only_sig = T,sep = "_")


rot_df <- as.data.frame(pca$rotation)
rot_df$label = rownames(rot_df)
rot_df %<>% tidyr::separate(label,into = c("ctx","type"),sep = "_")

p0 = rot_df %>% tidyr::gather("Sig","Weight",-ctx,-type) %>%
  dplyr::filter(Sig %in% glue::glue("PC{1:5}")) %>%
  ggplot(aes(color = ctx,
             x = forcats::fct_relevel(ctx,CMCP::pos_ms96),
             y=Weight)) +
  facet_grid(Sig~type) +
  geom_bar(stat = "identity",fill = NA) +
  scale_color_manual(values = CMCP::colors96) +
  theme(legend.position = "none")

p1 = rot_df %>%
  ggplot(aes(x = PC1,y=PC2,shape = type,color = ctx)) +
  geom_point() + scale_color_manual(values = CMCP::colors96) +
  geom_path(aes(group=ctx)) +
  theme(legend.position = "none")

samples_df = as.data.frame(pca$x)
samples_df$sample = rownames(pca$x)
samples_df$disease = pc[as.character(samples_df$sample)] %>% unlist()
samples_df$burden = bd[as.character(samples_df$sample)] %>% unlist()

p2 = samples_df %>% tidyr::gather("Sig","Weight",-sample,-disease,-burden) %>%
  dplyr::filter(Sig %in% glue::glue("PC{1:5}")) %>%
  ggplot(aes(x = disease, y = Weight, color = log10(burden))) +
  #geom_jitter() +
  geom_boxplot(fill = NA) +
  facet_grid(Sig~.,scales = "free") +
  theme(axis.text.x = element_text(angle = 90))


p3 = samples_df %>% ggplot(aes(x =PC1, y = PC2, color = disease)) +
  geom_point()

fp = plot_grid(p3,p1,p2,p0,labels = "AUTO")
