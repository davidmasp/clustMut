library(magrittr)

ips_tcia = readr::read_tsv("~/data/TCIA_IPS/patientsAll.tsv")
sample_sig = readRDS("~/projects/CMC/clustMut/output/w.rds")

sample_sig$barcode = sample_sig$sample %>%
  stringr::str_split("-") %>%
  purrr::map(function(x) paste0(x[1:3],collapse = "-")) %>%
  unlist()

matched_samples = sample_sig$barcode %in% ips_tcia$barcode %>% sum()
matched_per = scales::percent(matched_samples/nrow(sample_sig))
print(glue::glue("Matched samples in TCIA - {matched_per}"))


dat = ips_tcia %>% dplyr::inner_join(sample_sig) %>%
  dplyr::select(barcode,
                disease,
                burden,
                dplyr::starts_with("ips"),
                dplyr::starts_with("V"))

plot_list = list()

for (i in c("ips_ctla4_pos_pd1_pos",
            "ips_ctla4_neg_pd1_pos",
            "ips_ctla4_pos_pd1_neg",
            "ips_ctla4_neg_pd1_neg")){
  dat$ips_categories = plyr::mapvalues(as.data.frame(dat)[,i],
                                       from=c(3,4,5,6, 7, 8,9,10),
                                       to=c("low","low","low", "low", "medium","medium","high","high"))

  dat$ips_categories = factor(dat$ips_categories,levels = c("low","medium","high"))

  plot_list[[i]] = clustMut::multi_values_plot(dat = as.data.frame(dat),
                              category = "ips_categories",
                              values = c(glue::glue("V{1:8}"),"burden"))
}



