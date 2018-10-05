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
    type = "integer",
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
deps = c("genomicHelpersDMP",
         "magrittr",
         "readr",
         "purrr",
         "ggplot2",
         "viridis",
         "VariantAnnotation",
         "parallel",
         "glue",
         "clustMut")

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
  opt$data = "tests_exec/data/"
  opt$recursive = TRUE
  opt$glob = "*.tsv"

  opt$keepMSM = F
  opt$mutlist = F
  opt$keepVR = T
  opt$unclustkeep = FALSE

  opt$dist_cutoff = 500
  opt$fdr_method = "fdr"
  opt$nmuts = 1
  opt$outuput_prefix = glue::glue("clust_{opt$nmuts}")

  opt$verbose = T
}


print(opt)

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


# DISTANCE ===============================================================

# first we iterate over each file to compute the FDR, this means
# that it is not possible to handle a sample which is divided
# in 2 files. (CAVEAT)
vr_res = lapply(file_paths,function(x){
  #browser()
  dat = suppressMessages(
    readr::read_tsv(x,
                    col_types = list(
                      chr = readr::col_character()
                      )
                    )
    )
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
                      method = opt$fdr_method,
                      dist_cutoff = opt$dist_cutoff,
                      n = opt$nmuts)

  return(vr_res)
})
names(vr_res) = NULL
vr_res = do.call("c",vr_res)

if(opt$fdr_method == "fdr"){
  print("ploting files.")
  plot_exp(vr_res,filename = glue::glue("{opt$outuput_prefix}_plot.pdf"))
} else if (opt$fdr_method == "FDR"){
  print("ploting for FDR not implemented yet")
}

if (opt$verbose){
  all_samples = length(unique(VariantAnnotation::sampleNames(vr_res)))
  cat("\n")
  cat(cli::boxx(glue::glue("Number of samples detected: {all_samples}")))
  cat("\n")
}

# EVENTS ====================================================================

# we sort to be sure this will make only adjacent groups
vr_res = sortSeqlevels(vr_res)
vr_res = sort(vr_res)
if (opt$fdr_method == "fdr"){

  # we set up the recieving column
  vr_res$event_type = NA

  # we split again by sample and chr
  vr_res_list = split(vr_res,list(sampleNames(vr_res),seqnames(vr_res)))

  # we search events
  vr_res_list = lapply(vr_res_list, function(vr){
    vr$event_type = detect_events(vr$fdr,opt$fdr_cutoff,5,"kataegis")
    return(vr)
  })

  vr_res = unlist_GR_base_list(vr_res_list)

  omikli_mask = vr_res$fdr < opt$fdr_cutoff & is.na(vr_res$event_type)
  if (sum(omikli_mask) != 0){
    vr_res[omikli_mask]$event_type = "omikli"
  }

}



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

