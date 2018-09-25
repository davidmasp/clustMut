

stopifnot(requireNamespace("optparse",quietly = T))
stopifnot(requireNamespace("cli",quietly = T))
stopifnot(requireNamespace("clustMut",quietly = T))


cat(cli::rule(center = "ClustMut"))
cat("\n")
cat(cli::boxx(c(glue::glue("clustmut version: {packageVersion('clustMut')}"),
            "DMP: david.mas@irbbarcelona.org"),
          padding = 1, align = "center"))
cat("\n")
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
    c("-t", "--true_positive"),
    action = "store_true",
    default = FALSE,
    help = "tag to extrapolate spectrum to total number of true positives"
  ),
  make_option(
    c("-f", "--fdr_cutoff"),
    action = "store",
    default = 0.2,
    type = 'numeric',
    help = "fdr cutoff to compute the mutation subtype spectrum [default %default]"
  ),
  make_option(
    c("-n", "--cores"),
    action = "store",
    default =  NULL,
    type = 'integer',
    help = "Number of cores to parallelize [default %default]"
  ),make_option(
    c("-p", "--pairs_size"),
    action = "store",
    default = 10000,
    type = 'integer',
    help = "(vaf mode) Window size (in bp) used to find mutation pairs [default %default]"
  ),
  make_option(
    c("-S", "--simulation_size_input"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "(vaf mode) Number of max simulations of VAF [default %default]"
  ),
  make_option(
    c("-r", "--recursive"),
    action = "store_true",
    dest = "recursive",default = FALSE,
    help = "Make the program look inside folders within the data folder "
  ),
  make_option(
    c("-u", "--keep_uncl"),
    action = "store_true",
    dest = "unclustkeep",default = FALSE,
    help = "Keep also a unclustered MSM in the object file "
  ),
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    dest = "verbose",default = FALSE,
    help = "Give a more vebose output"
  ),
  make_option(
    c("-o", "--outuput_prefix"),
    action = "store",
    default = "output",
    type = 'character',
    help = "Prefix to base the output naming [default %default]"
  ),
  make_option(
    c("-M", "--mode"),
    action = "store",
    default = "distance",
    type = 'character',
    help = "Distance mode, to use the distance as statistic [default %default]"
  ),
  make_option(
    c("-N", "--nmuts"),
    action = "store",
    default = 1,
    help = "Number of mutations to consider as pairs [default %default]"
  ),
  make_option(
    c("-m", "--fdr_method"),
    action = "store",
    default = "fdr",
    type = 'character',
    help = "Fdr mode, fdr for local FDR for tail-based [default %default]"
  ),
  make_option(
    c("-d", "--dist_cutoff"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Distance cutoff for tail-based FDR [default %default]"
  ),
  make_option(
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
    c("-R", "--reference_genome"),
    action = "store",
    default = "Hsapiens.UCSC.hg19",
    help = "The genome reference used to compute MSM"
  )
)


opt = parse_args(OptionParser(option_list=option_list))

if (opt$verbose) {
  message(print(opt))
}

# check dependencies ===========================================================
stopifnot(requireNamespace("clustMut",quietly = T))
stopifnot(requireNamespace("genomicHelpersDMP",quietly = T))
stopifnot(requireNamespace("magrittr",quietly = T))
stopifnot(requireNamespace("readr",quietly = T))
stopifnot(requireNamespace("purrr",quietly = T))
stopifnot(requireNamespace("ggplot2",quietly = T))
stopifnot(requireNamespace("VariantAnnotation",quietly = T))


if (opt$verbose){
  library(VariantAnnotation)
  library(clustMut)
  library(genomicHelpersDMP)
  library(magrittr)
  if (!is.null(opt$cores)){
    stopifnot(requireNamespace("parallel",quietly = T))
    library(parallel)
  }
} else {
  suppressPackageStartupMessages(library(VariantAnnotation))
  suppressPackageStartupMessages(library(clustMut))
  suppressPackageStartupMessages(library(genomicHelpersDMP))
  suppressPackageStartupMessages(library(magrittr))
  if (!is.null(opt$cores)){
    stopifnot(requireNamespace("parallel",quietly = T))
    suppressPackageStartupMessages(library(parallel))
  }
}

# reading files ================================================================

if (interactive()){
  opt$mode = "distance"
  opt$data = "Y:/users/dmas/data/ICGC_MUTS/MELA_AU/splited/"
  opt$recursive = TRUE
  opt$glob = "*-randomized.tsv"

  opt$keepMSM = T
  opt$mutlist = T
  opt$keepVR = T
  opt$unclustkeep = TRUE

  opt$dist_cutoff = 500
  opt$fdr_method = "FDR"
  opt$nmuts = 1
  opt$outuput_prefix = glue::glue("clust_{opt$nmuts}")

  opt$verbose = T
  opt$true_positive = T

  opt$reference_genome = "Scerevisiae.NIEHS.ySR127"
}

# opt$mode = "vaf"
# opt$data = "E:/local-data/TCGA_MUTS/TCGA_VR/"
# opt$glob = "*_VR.rds"
# opt$recursive = TRUE
#  opt$alignability_mask = "E:/local-data/CRG_alignability/hg19/LEGACY/crg36AlignExtNoBlackRmvsh19_RngMask_savedInt=TRUE.bed"

# opt$mode = "edit"
# opt$data = "E:/local-data/TCGA_MUTS/TCGA_VR/"
# opt$glob = "*_VR.rds"
# opt$recursive = TRUE
#  opt$alignability_mask = "E:/local-data/CRG_alignability/hg19/LEGACY/crg36AlignExtNoBlackRmvsh19_RngMask_savedInt=TRUE.bed"

# MAIN ========================

path = opt$data
file_paths = fs::dir_ls(path,
                        glob = opt$glob,
                        recursive = opt$recursive)

if (length(file_paths) == 0){
  stop("No files provided.")
}

if (opt$mode == "distance"){ # if distance -i should be randommut out
  # DISTANCE ===============================================================
  library(progress)
  pb <- progress_bar$new(
    format = "Distance analysis, :percent done. eta: :eta",
    total = length(file_paths), clear = FALSE)

  vr_res = lapply(file_paths,function(x){

    dat = suppressMessages(readr::read_tsv(x))
    original = nrow(dat)

    # apparently this happens in ICGC
    mask = (dat$ref == dat$alt)
    dat = dat[!mask,]

    # check for duplicates
    dat %<>% dplyr::distinct(chr,start,sample,alt,.keep_all=TRUE)

    # remove NNN and outside set
    ref_in = genomicHelpersDMP::dna_codes[[opt$pair_set]]
    mask=dat$ref %in% ref_in

    dat = dat[mask,]

    after_filter = nrow(dat)
    per = (1 - (after_filter/original)) %>% scales::percent()
    if (opt$verbose){
      print(glue::glue("{per} mutations discarded in sample {x}."))
    }

    # analysis ===================
    tmp = parse_randommut_vr(dat)

    vr_res = clust_dist(vr = tmp$VR,
                        rand_df = tmp$RAND,
                        no_cores = opt$cores,
                        method = opt$fdr_method,
                        dist_cutoff = opt$dist_cutoff,
                        n = opt$nmuts)
    pb$tick()
    return(vr_res)
  })
  names(vr_res) = NULL
  vr_res = do.call("c",vr_res)

  if(opt$fdr_method == "fdr"){
    print("ploting files.")
    plot_exp(vr_res,filename = glue::glue("{opt$outuput_prefix}_plot.pdf"))
  }

  all_samples = length(unique(VariantAnnotation::sampleNames(vr_res)))
  cat("\n")
  cat(cli::boxx(glue::glue("Number of samples detected: {all_samples}")))
  cat("\n")

} else if (opt$mode == "vaf") { # read rds VR files
  # VAF ====================================================================
  dat = dat = VR_preprocessing(file_paths = file_paths,
                               pair_set = opt$pair_set,
                               alignability_mask = opt$alignability_mask
  )

  # analysis ===================
  # This could be parallelized
  pb <- progress_bar$new(
    format = "Computing VAF fdr :percent eta: :eta",
    total = length(file_paths), clear = FALSE)


  vr_res = purrr::map(dat,function(vr){
    vr_res = clustMut::local_vaflr_fdr(vr = vr,
                                       simulation_size_input =
                                       opt$simulation_size_input,
                                       pairs_size = opt$pairs_size)

    pb$tick()
    return(vr_res)
  })

  names(vr_res) = NULL
  vr_res <- do.call(what = "c", args = vr_res)
  all_samples = length(unique(VariantAnnotation::sampleNames(vr_res)))
  cli::boxx(glue::glue("Number of samples detected: {all_samples}"))




} else if (opt$mode == "edit"){
  # EDIT ===================================================================

  # read rds VR file
  dat = VR_preprocessing(file_paths = file_paths,
                         pair_set = opt$pair_set,
                         alignability_mask = opt$alignability_mask
                         )

  # # analysis EDIT ===================
  # # This could be parallelized
  # pb <- progress_bar$new(
  #   format = "Computing Edit distance fdr :percent eta: :eta",
  #   total = length(file_paths), clear = FALSE)


  # Initiate cluster
  print("Running a parallelized version of the edit dist.")
  library(parallel)
  cl <- makeCluster(opt$cores)
  clusterEvalQ(cl = cl,library("clustMut"))
  clusterEvalQ(cl = cl,library("genomicHelpersDMP"))
  clusterExport(cl = cl,varlist = c("opt"),envir=environment())

  vr_res = parLapply(cl = cl,X = dat,fun = function(vr){

    genome = genome_selector()
    vr_res = edit_distance_fdr(vr = vr,
                               k = 4,
                               genome = genome,
                               pairs_size = opt$pairs_size,
                               simulation_size = opt$simulation_size_input )
    #pb$tick()
    return(vr_res)
  })

  stopCluster(cl)

  names(vr_res) = NULL
  vr_res <- do.call(what = "c", args = vr_res)
  all_samples = length(unique(VariantAnnotation::sampleNames(vr_res)))
  cli::boxx(glue::glue("Number of samples detected: {all_samples}"))
} else {
  stop("Method not implemented")
}


# output ==============================================================


return_output(
  vr_res = vr_res,
  keepVR = opt$keepVR,
  outuput_prefix = opt$outuput_prefix,
  mode = opt$mode,
  fdr_cutoff = opt$fdr_cutoff,
  mutlist = opt$mutlist,
  keepMSM = opt$keepMSM,
  unclustkeep = opt$unclustkeep,
  kmer = opt$kmer,
  true_positive = opt$true_positive,
  reference_genome = opt$reference_genome,
  fdr_method = opt$fdr_method,
  dist_cutoff = opt$dist_cutoff
)

