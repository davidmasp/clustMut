# CLUSTMUT

stopifnot(requireNamespace("optparse",quietly = T))
stopifnot(requireNamespace("cli",quietly = T))
stopifnot(requireNamespace("clustMut",quietly = T))

# optparse =====================================================================
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(
    c("-Y","--yawn"),
    action = "store_true",
    default = FALSE,
    help = "Print version of the package and quit",
    dest = "version"
  ),
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
deps = c("clustMut",
         "genomicHelpersDMP",
         "magrittr",
         "readr",
         "purrr",
         "ggplot2",
         "VariantAnnotation",
         "glue",
         "parallel")

stopifnot(requireNamespace("clustMut",quietly = T))
clustMut::clustmut_internal_check_dependencies(deps)

# print version ===============================================================
clustMut::clustmut_internal_return_version(version = (opt$verbose | opt$version),
               quit = opt$version)

# load packages ===============================================================
libs = c(
  "VariantAnnotation",
  "clustMut",
  "genomicHelpersDMP",
  "magrittr",
  "parallel"
)

stopifnot(all(libs %in% deps))

if (opt$verbose){
  for (i in libs){
    library(i,character.only = TRUE)
  }
} else {
  for (i in libs){
    suppressPackageStartupMessages(library(i,character.only = TRUE))
  }
}

# interactive debugging options =============================================
if (interactive()){
  opt$mode = "distance"
  opt$data = "tests_exec/data/"
  opt$recursive = TRUE
  opt$glob = "*.tsv"

  opt$keepMSM = F
  opt$mutlist = F
  opt$keepVR = T
  opt$unclustkeep = TRUE

  opt$dist_cutoff = 500
  opt$fdr_method = "FDR"
  opt$nmuts = 1
  opt$outuput_prefix = glue::glue("clust_{opt$nmuts}")

  opt$verbose = T
  opt$reference_genome = "Scerevisiae.NIEHS.ySR127"
}

# MAIN ========================

# reading files ================================================================
path = opt$data
file_paths = fs::dir_ls(path,
                        glob = opt$glob,
                        recursive = opt$recursive)

# check if files are available and exists
if (length(file_paths) == 0){
  stop("No files provided.")
}


# VAF ====================================================================
dat = VR_preprocessing(file_paths = file_paths,
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


# output ==============================================================


clustmut_internal_return_output(
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

