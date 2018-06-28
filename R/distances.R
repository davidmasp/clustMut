# =========================================================================== #
#                      _____  _     _
#                     |  __ \(_)   | |
#                     | |  | |_ ___| |_ __ _ _ __   ___ ___
#                     | |  | | / __| __/ _` | '_ \ / __/ _ \
#                     | |__| | \__ \ || (_| | | | | (_|  __/
#                     |_____/|_|___/\__\__,_|_| |_|\___\___|
#
# ==================================================================+======== #



# ======= LOW LEVEL ======= #

#' Get mutation distance
#'
#' Get the distance between a particular mutaion and the next or to the previous.
#' @param x
#' @param k
#'
#' @return
#' @export
#'
#' @examples
compute_distances <- function(x,k = 1,use = dplyr::lag){
  original_order <- order(x)
  x = x[original_order]
  distance <- x - use(x,n=k)
  distance <- abs(distance)

  # OPEN A ISSUE - need to check if this breaks smth
  # if(all(is.na(distance[1:k]))){
  #   # open a issue for this
  #   distance[1:k] <- x[1:k]
  # } else if (all(is.na(distance[(length(x)-k):length(x)]))){
  #   distance[(length(x)-k):length(x)-1] = x[length(x)-1] - x[length(x)]
  # }

  distance <- distance[order(original_order)]
  return(distance)
}

#' Get minimum adjacent distance
#'
#' Wrapper of compute_distances
#' Get the distance between each mutation and the closest adjacent one.
#' @param x vector of positions (doesn't check same chromosome)
#' @param k adjacent positions (k = 1 by default)
#' @param use function to select appropiate adjacent mutation (min by default)
#'
#' @return
#' @export
#'
#' @examples
compute_m_distance <- function(x,k=1,use = min){
  prev = compute_distances(x,k=k,use = dplyr::lag)
  post = compute_distances(x,k=k,use = dplyr::lead)

  # this is to remove NA
  prev[is.na(prev)] = post[is.na(prev)]
  post[is.na(post)] = prev[is.na(post)]

  # to check this matrix(c(c(1,2,3,2),c(1,1,1,1)),ncol = 2)
  dist = matrix(c(prev,post),ncol = 2)
  m_dist = apply(dist,1,use)
  return(m_dist)
}


# ======= PER SAMPLE LEVEL ======= #

#' Randomized distances
#'
#' Compute minimum distance at each column of the data frame.
#'
#' @param x a dataframe or tibble object
#' @param f factor to group rows, if more than one group is required you should input a list containing both factors.
#' @param k number of mutations that should be enclosed in the mutation pair.
#' @param no_cores if null normal execution, if integer, parallel mode is triggered with the number of cores specified.
#'
#' @return
#' @export
#'
#' @examples
compute_distances_splited_tbl <- function(x,
                                          f,
                                          k=1,
                                          no_cores = NULL){
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("purrr", quietly = TRUE)
  # prepare cluster
  if (!is.null(no_cores)){
    #browser()
    requireNamespace("parallel", quietly = TRUE)
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterEvalQ(cl, { library(clustMut)})
  }


  if (is.list(f)) f <- interaction(f,drop = TRUE)

  if (!is.null(no_cores)){
    df_splited = base::split(x = x,f = f)
    rand_dist = parallel::parLapply(cl = cl,
                                    X = df_splited,
                                    fun = function(x){
        dist_df = apply(x,2,compute_m_distance,k = k)
        if(nrow(x) == 1){
          warning("1 position in group")
          dist_df = t(dist_df)
        }
        dist_df = data.frame(dist_df)
        return(dist_df)
      })
  } else {

    rand_dist = x %>% split(f) %>% purrr::map(function(x){

      dist_df = apply(x,2,compute_m_distance,k=k)

      if(nrow(x) == 1){
        warning("1 position in group")
        dist_df = t(dist_df)
      }

      dist_df = data.frame(dist_df)
      return(dist_df)
    })
  }

  rand_dist = dplyr::bind_rows(rand_dist)

  if (!is.null(no_cores)){
    parallel::stopCluster(cl)
  }

  return(rand_dist)
}




#' Find IMD per mutation k pairs
#'
#' In a VR object, finds mutation pairs that enclose n number of mutations. It
#' only selects the nearest pair.
#' As a general way to handle this I would use nearesDistance if k=1.
#'
#' @param vr a VRanges object
#' @param enclosing Number of mutations between pairs
#'
#' @return a VRanges object with a mdist column
#' @export
#'
#' @examples
compute_distance_vr <- function(vr,enclosing){
  stopifnot(length(unique(sampleNames(vr))) == 1)
  # single sample assumption
  ## optional, if you want a genomic order of the chromosomes
  stopifnot(!is.unsorted(vr))

  mcols(vr)$distance = NA
  ## split into a GRangesList
  ## where each element has all ranges for one chromosome
  vrl = split(vr, seqnames(vr))

  ## apply a function to the ranges of each chromosome
  res = lapply(names(vrl), function(x){
    VR=vrl[[x]]
    stopifnot(length(unique(seqnames(VR))) == 1)
    dist = compute_m_distance(x = start(VR),k = enclosing,use = min)
    mcols(VR)$distance = dist
    return(VR)
  })

  # merge the result back
  names(res) <- NULL
  res = do.call("c",res)
  return(res)
}



