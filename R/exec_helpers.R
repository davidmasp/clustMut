# Exec helpers ===================
# Functions that should go in the executable files.
# Output, pruning etc.

# THIS SHOULD NOT BE EXPORTED! ISSUE



#' Check dependencies
#'
#' @param deps character vector with package names
#'
#' @return
#' @export
#'
#' @examples
clustmut_internal_check_dependencies <- function(deps) {
  for (i in deps){
    if (!requireNamespace(i,quietly = T)){
      cat(paste("Dependency",i,"is missing in the installation"))
    }
  }
}


#' Return version of code
#'
#' @param version TRUE if you want to return the version of the code
#' @param quit TRUE if quitting after the version has been rinted in the CLI
#'
#' @return NULL
#' @export
#'
#' @examples
clustmut_internal_return_version <- function(version,quit){
  if (version){
    cat(cli::rule(center = "ClustMut"))
    cat("\n")
    cat(cli::boxx(c(glue::glue("clustmut version: {packageVersion('clustMut')}"),
                    "DMP: david.mas@irbbarcelona.org"),
                  padding = 1, align = "center"))
    cat("\n")
    if (quit){
      stop("Quitting...")
    }

  } else{
    return(NULL)
  }
}


clustmut_internal_return_output <- function(vr_res,
                          keepVR,
                          outuput_prefix,
                          mode,
                          fdr_cutoff,
                          mutlist,
                          keepMSM,
                          unclustkeep,
                          kmer,
                          true_positive,
                          reference_genome,
                          fdr_method,
                          dist_cutoff){
  # OUTPUT ==================================================================

  print("Print Output")

  # save the VR object

  if (keepVR){
    saveRDS(object = vr_res,
            file = glue::glue("{outuput_prefix}_{mode}_VRanges.rds"))
  }

  if (fdr_method == "fdr"){
    selected_muts = vr_res[vr_res$fdr < fdr_cutoff & !is.na(vr_res$fdr)]
  } else if (fdr_method == "FDR"){
    print(colnames(vr_res))
    warning("this has not yet been implemented, returning NULL")
    return(NULL)
    selected_muts = vr_res[vr_res$fdr < fdr_cutoff & !is.na(vr_res$fdr)]
  }



  # save a list output
  if (mutlist){

    mutid = paste(seqnames(selected_muts),
                  start(selected_muts),
                  sampleNames(selected_muts),
                  alt(selected_muts),
                  sep = ":")

    readr::write_lines(mutid,path = glue::glue("{outuput_prefix}_{mode}_mutlist.txt"))
  }


  if (keepMSM){
    # compute the MSM -> see issue #5
    # this could also be parallelized too.
    MSM_clust = compute_MSM(vr = selected_muts,
                            k = kmer,
                            tp = true_positive,
                            genome = genome_selector(
                              alias = reference_genome))


    if (unclustkeep){
      vr_unclust = vr_res[vr_res$fdr>=fdr_cutoff | is.na(vr_res$fdr) ]
      MSM_uncl = compute_MSM(vr = vr_unclust,
                             k = kmer,
                             tp = true_positive,
                             genome = genome_selector(
                               alias = reference_genome) )
      # I output all samples, because it is a list.
      MSM_result = list(clust = MSM_clust,uncl = MSM_uncl)
    } else {
      MSM_result = list(clust = MSM_clust)
    }

    saveRDS(object = MSM_result, file = glue::glue("{outuput_prefix}_{mode}_MSM.rds"))
  }
}
