###############################################################################
#                         _           _                
#                        | |         | |               
#                     ___| |_   _ ___| |_ ___ _ __ ___ 
#                    / __| | | | / __| __/ _ \ '__/ __|
#                   | (__| | |_| \__ \ ||  __/ |  \__ \
#                    \___|_|\__,_|___/\__\___|_|  |___/
#
###############################################################################
                                    
                                    
compute_fdr_parallel <- function(dat_gr,rand_dist,f,no_cores){
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("purrr", quietly = TRUE)
  requireNamespace("localFDR", quietly = TRUE)
  requireNamespace("broom", quietly = TRUE)
  #browser()
  df = data.frame(mdist = dat_gr$mdist)
  df = cbind(df,data.frame(rand_dist))
  dat_split = base::split(df,f = f)

  cl <- parallel::makeCluster(no_cores)
  parallel::clusterEvalQ(cl, {
    library(clustMut)
    library(magrittr)
    library(localFDR)
    })

  fdr_vec = parallel::parLapply(cl = cl,X = dat_split,fun = function(x){
    pos_distance = x$mdist
    random_matrix = x %>% dplyr::select(-mdist) %>% as.matrix()

    # this is to avoid log(0) should be done internally I guess.
    pos_distance = pos_distance + 1
    random_matrix = random_matrix + 1

    fdr_matrix = apply(random_matrix,2,function(y){
      localFDR::compute_densityfdr(obs = log(pos_distance), null = log(y))
    })

    fdr_vec = apply(fdr_matrix,1,median)

    return(fdr_vec)
  })

  stopCluster(cl)

  dat_gr$fdr =  unlist(fdr_vec)

  return(dat_gr)
}



compute_fdr_default <- function(dat_gr,rand_dist,f){
  requireNamespace("purrr", quietly = TRUE)
  requireNamespace("localFDR", quietly = TRUE)
  requireNamespace("broom", quietly = TRUE)
  pos_split = base::split(dat_gr$mdist,f = f)
  rand_split = split(rand_dist,f=f)
  #browser()
  lol = list(x = pos_split,y = rand_split,k = names(pos_split), z = names(rand_split))
  fdr_vec = purrr::pmap(.l = lol,.f = function(x,y,k,z){
    if(k != z){stop("lists do not match")}

    pos_distance = x
    random_matrix = as.matrix(y)

    # this is to avoid log(0) should be done internally I guess.
    pos_distance = pos_distance + 1
    random_matrix = random_matrix + 1

    fdr_matrix = apply(random_matrix,2,function(y){
      localFDR::compute_densityfdr(obs = log(pos_distance), null = log(y))
    })

    fdr_vec = apply(fdr_matrix,1,median)

    return(fdr_vec)

  } )

  dat_gr$fdr =  unlist(fdr_vec)

  return(dat_gr)
}

compute_fdr <- function(dat_gr,rand_dist,f,no_cores=NULL){

  if (is.list(f)) f <- interaction(f, drop = TRUE)

  if (!is.null(no_cores)){
    dat_fdr = compute_fdr_parallel(dat_gr,rand_dist,f,no_cores = no_cores)
  } else {
    dat_fdr = compute_fdr_default(dat_gr,rand_dist,f)
  }

  return(dat_fdr)

}
