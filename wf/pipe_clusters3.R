source("wf/snip/load_cmc.R")

a = genomicHelpersDMP::list_matrix_sample(matrix_list = MSM_list,
                                      reference_matrix = "clust")

library(NMF)
nmf_test = NMF::nmf(a,2:20,nrun = 50)
f_silouette = plot(nmf_test)

selected_rank = 8

nmf_fit = NMF::nmf(a,selected_rank,nrun = 200)

nmf_sigs = as.data.frame(t(nmf_fit@fit@H))
nmf_sigs$label = rownames(nmf_sigs)
nmf_sigs %<>% tidyr::separate(label,into = c("ctx","type"),sep = "_")

theme_set(theme_cowplot(font_size = 8))
pbase = nmf_sigs %>% tidyr::gather("Sig","Weight",-ctx,-type) %>%
  ggplot(aes(fill = ctx,
             x = forcats::fct_relevel(ctx,CMCP::pos_ms96),
             y=Weight)) +
  facet_grid(Sig~type,scales = "free",switch = "x") +
  geom_bar(stat = "identity",width = 1) +
  scale_fill_manual(values = CMCP::colors96) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90))

pleft_df = nmf_sigs %>% tidyr::gather("Sig","Weight",-ctx,-type) %>%
  tidyr::spread(key=type,value = Weight)
pleft_df$cte = "cte"
pleft = pleft_df %>%  ggplot(aes(y = clust,x = uncl,
             color = ctx)) +
  geom_point() +
  facet_wrap(Sig~cte,scales = "free",ncol = 1) +
  geom_abline(slope = 1,intercept = 0) +
  scale_color_manual(values = CMCP::colors96) +
  theme(legend.position = "none",
        strip.background  = element_blank(),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


fp = plot_grid(pleft,pbase,rel_widths = c(.5,9.5),align="h",axis = "b")



w_df = nmf_fit@fit@W %>% as.data.frame()
w_df$sample = rownames(w_df)
w_df %<>% tibble::as.tibble()
w_df$disease = pc[w_df$sample] %>% unlist()
w_df$burden = bd[w_df$sample] %>% unlist()
w_df %<>% tidyr::separate(disease,into = c("CT","Center"),sep="-")


burden_plot = w_df %>% tidyr::gather("Sig","Weight",-sample,-CT,-Center,-burden) %>%
  ggplot(aes(x = burden,y=Weight,color = CT)) +
  geom_point(alpha = 0.3) +
  scale_y_sqrt() +
  scale_x_log10() +
  facet_grid(Sig~CT,scales = "free") +
  geom_smooth(method = "lm" ) +
  theme(legend.position = "none")


top_row = plot_grid(f_silouette+theme_cowplot(font_size = 8)+theme(legend.position = "bottom"),burden_plot,rel_widths = c(4,6),labels = c("A","B"))
plot_grid(top_row,fp,nrow = 2,labels = c(NA,"C"))
ggsave("test.pdf",height = 12,width = 17)
