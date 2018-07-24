
stopifnot(requireNamespace("optparse",quietly = T))
stopifnot(requireNamespace("cli",quietly = T))
stopifnot(requireNamespace("clustMut",quietly = T))


cli::rule(center = "ClustMut - Custom method")
cli::boxx(c(glue::glue("clustmut version: {packageVersion('clustMut')}"),
            "DMP: david.mas@irbbarcelona.org"),
          padding = 1, align = "center")



# optparse =====================================================================
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-i", "--data"),
    action = "store",
    type = 'character',
    help = "Data folder"
  ),
  make_option(
    c("-r", "--recursive"),
    action = "store_true",
    dest = "recursive",default = FALSE,
    help = "Make the program look inside folders within the data folder "
  ),
  make_option(
    c("-a", "--alignability_mask"),
    action = "store",
    type = 'character',
    help = "(vaf mode) Bed file containing the included regions. Mutations there will be included [default %default]",
    default = NULL
  ),
  make_option(
    c("-s", "--pair_set"),
    action = "store",
    default = "N",
    type = 'character',
    help = "Mutation set to filter mutations before calling clusters [default %default]"
  ),
  make_option(
    c("-g", "--glob"),
    action = "store",
    default = "*-randomized.tsv",
    type = 'character',
    help = "GLOB used to look for each individual file. [default %default]"
  ),
  make_option(
    c("-k", "--kmer"),
    action = "store",
    default = 1,
    type = 'integer',
    help = "Number of downstream nucleotides to extract the MS [default %default]"
  ),
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    dest = "verbose",default = FALSE,
    help = "Give a more vebose output"
  ),
  make_option(
    c("-I", "--IMD"),
    action = "store",
    dest = "imd",default = 10000,
    help = "Cluster length"
  ),
  make_option(
    c("-N", "--nmuts"),
    action = "store",
    dest = "nmuts",default = 6,
    help = "Number of mutations in a cluster"
  ),
  make_option(
    c("-n", "--cores"),
    action = "store",
    default = 5,
    type = 'integer',
    help = "Number of cores to parallelize [default %default]"
  ),make_option(
    c("-V", "--keepVR"),
    action = "store_true",
    default = FALSE,
    help = "Keep the resulting VRanges object with the FDR/fdr values [default %default]"
  ),
  make_option(
    c("-l", "--mutlist"),
    action = "store_true",
    default = FALSE,
    help = "Output the selected mutations as a list. [default %default]"
  ),
  make_option(
    c("-w", "--keepMSM"),
    action = "store_true",
    default = FALSE,
    help = "Compute and output the mutation subtype matrix. [default %default]"
  ),
  make_option(
    c("-o", "--outuput_prefix"),
    action = "store",
    default = "output",
    type = 'character',
    help = "Prefix to base the output naming [default %default]"
  ),
  make_option(
    c("-u", "--keep_uncl"),
    action = "store_true",
    dest = "unclustkeep",default = FALSE,
    help = "Keep also a unclustered MSM in the object file "
  ),
  make_option(
    c("-R", "--reference_genome"),
    action = "store",
    default = "Hsapiens.UCSC.hg19",
    help = "The genome reference used to compute MSM"
  )
)



## reading



opt = parse_args(OptionParser(option_list=option_list))

if (opt$verbose) {
  message(print(opt))
}

# check dependencies ===========================================================
stopifnot(requireNamespace("clustMut",quietly = T))
stopifnot(requireNamespace("VariantAnnotation",quietly = T))

if (opt$verbose){
  library(clustMut)
  library(genomicHelpersDMP)
  library(magrittr)
  if (!is.null(opt$cores)){
    stopifnot(requireNamespace("parallel",quietly = T))
    library(parallel)
  }
} else {
  library(clustMut)
  library(genomicHelpersDMP)
  library(magrittr)
  if (!is.null(opt$cores)){
    stopifnot(requireNamespace("parallel",quietly = T))
  }
}


if (interactive()){
  opt$data= "~/data/TCGA_MUTS/TCGA_VR/LGG/"
  opt$glob = "*_VR.rds"
  opt$recursive = TRUE
  opt$alignability_mask = "~/data/CRG_alignability/hg19/LEGACY/crg36AlignExtNoBlackRmvsh19_RngMask_savedInt=TRUE.bed"

  opt$keepMSM = TRUE
  opt$keepVR = TRUE
  opt$mutlist = TRUE
  opt$keep_uncl = TRUE
}


path = opt$data
file_paths = fs::dir_ls(path,
                        glob = opt$glob,
                        recursive = opt$recursive)

dat = VR_preprocessing(file_paths = file_paths,
                       pair_set = opt$pair_set,
                       alignability_mask = opt$alignability_mask)




if (opt$cores == 1 | is.null(opt$cores)){
  library(progress)
  pb <- progress_bar$new(
    format = "Computing clusters by Custom method :percent eta: :eta",
    total = length(file_paths), clear = FALSE)

  vr_res = purrr::map(dat,function(vr){
    vr_res = clustMut::custom_basic_clustering(vr = vr,
                                               IMD = opt$imd,
                                               nmuts = opt$nmuts)

    pb$tick()
    return(vr_res)
  })

} else {
  cat("entering parallel mode")
  cl = parallel::makeCluster(opt$cores)

  parallel::clusterEvalQ(cl = cl,library(clustMut))
  clusterEvalQ(cl = cl,library(genomicHelpersDMP))
  clusterExport(cl = cl,varlist = c("opt"),envir=environment())

  vr_res = parLapply(cl = cl,
                     X = dat,
                     fun = function(vr){
                       vr_res = clustMut::custom_basic_clustering(vr = vr,
                                                                  IMD = opt$imd,
                                                                  nmuts = opt$nmuts)
                       return(vr_res)
                     })

  stopCluster(cl)
}



names(vr_res) = NULL
vr_res <- do.call(what = "c", args = vr_res)
all_samples = length(unique(VariantAnnotation::sampleNames(vr_res)))
cli::boxx(glue::glue("Number of samples detected: {all_samples}"))



# OUTPUT =======================================================================

print("Print Output")

# save the VR object


if (opt$keepVR){
  saveRDS(object = vr_res,
          file = glue::glue("{opt$outuput_prefix}_custom_w{opt$imd}_n{opt$nmuts}_VRanges.rds"))
}

selected_muts = vr_res[vr_res$custom_clust]


# save a list output
if (opt$mutlist){

  mutid = paste(seqnames(selected_muts),
                start(selected_muts),
                sampleNames(selected_muts),
                alt(selected_muts),
                sep = ":")

  readr::write_lines(mutid,path = glue::glue("{opt$outuput_prefix}_custom_w{opt$imd}_n{opt$nmuts}_mutlist.txt"))
}


if (opt$keepMSM){
  # compute the MSM -> see issue #5
  # this could also be parallelized too.

  MSM_clust = compute_MSM(vr = selected_muts,
                          k = opt$kmer,
                          tp = FALSE,
                          genome = genome_selector(
                            alias = opt$reference_genome))


  if (opt$unclustkeep){
    MSM_uncl = compute_MSM(vr = vr_res[!vr_res$custom_clust],
                           k = opt$kmer,
                           tp = FALSE,
                           genome = genome_selector(
                             alias = opt$reference_genome))

    MSM_uncl = MSM_uncl[rownames(MSM_clust),]


    MSM_result = list(clust = MSM_clust,uncl = MSM_uncl)
  } else {
    MSM_result = list(clust = MSM_clust)
  }

  saveRDS(object = MSM_result,
          file = glue::glue("{opt$outuput_prefix}_custom_w{opt$imd}_n{opt$nmuts}_MSM.rds"))
}




