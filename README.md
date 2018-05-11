# clustMut

The goal of clustMut is to ...

## Example

This is a basic example which shows you how to solve a common problem:

``` r

library(clustMut)

# steps single sample
path = Sys.getenv("path_to_sample")
dat = readr::read_tsv(path)
tmp = parse_randommut_vr(dat)
vr_res = clust_dist_sample(vr = tmp$VR,rand_df = tmp$RAND)
vr_res 
genomicHelpersDMP::compute_MSR(vr = vr_res[vr_res$fdr<0.2],
                               tp = unique(vr_res$tp))


# steps multiple
file_paths= fs::dir_ls(Sys.getenv("path_to_randomut_out"),
                       glob = "*_WGS_ssm_tcga_conf.tsv-randomized.tsv",
                       recursive = TRUE)
dat = purrr::map_df(file_paths,readr::read_tsv)
tmp = parse_randommut_vr(dat)
unique(VariantAnnotation::sampleNames(tmp$VR))
vr_res = clust_dist(vr = tmp$VR,rand_df = tmp$RAND,no_cores = 5)
vr_res
MSM = genomicHelpersDMP::compute_MSM(vr = vr_res[vr_res$fdr<0.2],
                               tp = T)


```
