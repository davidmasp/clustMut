# clustMut

The  goal of clustMut is to ...
 
## Example

This is a basic example which shows you how to solve a common problem:

``` r

library(clustMut)

path = Sys.getenv("file_path")
dat = readr::read_tsv(path)
dat = clustMut::parse_randommut_out(dat)
dat_gr = dat$GR
rand_df = dat$RAND

dat_gr_cl = sample_free_clusters(dat_gr,rand_df,plot = TRUE)
ggsave("test.pdf")
```
