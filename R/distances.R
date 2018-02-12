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
  distance <- x - use(x,k=k)
  distance <- abs(distance)
  if(is.na(distance[1])){
    distance[1] <- x[1]
  } else if (is.na(distance[length(x)])){
    distance[length(x)] = distance[length(x)-1]
  }

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
  prev = compute_distances(x,k=1,use = dplyr::lag)
  post = compute_distances(x,k=1,use = dplyr::lead)
  # to check this matrix(c(c(1,2,3,2),c(1,1,1,1)),ncol = 2)
  dist = matrix(c(prev,post),ncol = 2)
  m_dist = apply(dist,1,use)
  return(m_dist)
}


# ======= PER SAMPLE LEVEL ======= #

# order mantained in parallel execution
# I checked a bit and seems to be the case

#'Compute Distances in splited GR
#'
#' Compute minimum distance between adjacent mutation spliting by a factor. This
#' factor can be either the seqnames (chromosomes) and/or samples.
#'
#' @param gr Genomic Ranges object.
#' @param f factor to group rows, if more than one group is required you should input a list containing both factors.
#' @param no_cores if null normal execution, if integer, parallel mode is triggered with the number of cores specified.
#'
#' @return
#' @export
#'
#' @examples
compute_distances_splited <- function(gr,
                                     f,
                                     no_cores = NULL
                                     ){

  # prepare cluster
  if (!is.null(no_cores)){
    #browser()
    requireNamespace("parallel", quietly = TRUE)
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterEvalQ(cl, { library(clustMut)})
  }

  total_muts = length(gr)

  if (is.list(f)) f <- interaction(f, drop = TRUE)

  gr_per_chr =  GenomicRanges::split(gr,f)

  if (!is.null(no_cores)){
    mdist = end(gr_per_chr) %>%
      parallel::parLapply(cl = cl,X = .,fun = compute_m_distance) %>% unlist()
    gr$mdist = mdist
  } else {
    #browser()
    mdist =  unlist(lapply(X = end(gr_per_chr),FUN = compute_m_distance))
    gr$mdist = mdist
  }

  if (!is.null(no_cores)){
    parallel::stopCluster(cl)
  }

  return(gr)
}

#' Randomized distances
#'
#' Compute minimum distance at each column of the data frame.
#'
#' @param x a dataframe or tibble object
#' @param f factor to group rows, if more than one group is required you should input a list containing both factors.
#' @param no_cores if null normal execution, if integer, parallel mode is triggered with the number of cores specified.
#'
#' @return
#' @export
#'
#' @examples
compute_distances_splited_tbl <- function(x,
                                          f,
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
        dist_df = apply(x,2,compute_m_distance)
        if(nrow(x) == 1){
          warning("1 position in group")
          dist_df = t(dist_df)
        }
        dist_df = data.frame(dist_df)
        return(dist_df)
      })
  } else {

    rand_dist = x %>% split(f) %>% purrr::map(function(x){

      dist_df = apply(x,2,compute_m_distance)

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




