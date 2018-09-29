### fdr code


##### Density fdr ####

#' Transform a vector to monotonic
#'
#' @param x the vector
#' @param max maximum value allowed
#' @param ascending TRUE if vector should go from 0 to max, False if max to 0
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
to_monotonic <- function(x,
                         max = 1,
                         ascending = TRUE,
                         ...) {
  if (ascending) {
    res = pmin(max,cummax(x))
  } else if (!ascending) {
    res = pmin(max,cummin(x))
  } else{
    stop("jsdfkjsdhk")
  }
  return(res)
}


#' Gaussian density local fdr estimator
#'
#' It computes the local fdr estimate based on the density estimates from a nonnull and a null distribution.
#'
#' @param obs nonnull distribution of values
#' @param null expected or null distribution
#' @param monotonic if fdr should be monotonic
#' @param alternative left if fdr should be ascending, right if decending
#'
#' @return local fdr vector
#' @export
#'
#' @examples
compute_densityfdr <- function(obs,
                               null,
                               monotonic=TRUE,
                               alternative = c("left","right")){
  stopifnot(requireNamespace("broom",quietly = TRUE))
  require(magrittr)

  # i compute the extremes of both distributions
  ext = range(c(obs,null))

  # then I compute the densities using density() function
  obs_df = density(obs,from = ext[1],to=ext[2]) %>% broom::tidy()
  obs_df$type = "obs"

  exp_df = density(null,from = ext[1],to=ext[2]) %>% broom::tidy()
  exp_df$type = "null"

  # i check the x values are the same
  if (!all(obs_df$x == exp_df$x)){stop("Jjhdskjbfkjds")}

  # i compute the local fdr per density estimated point
  fdr = data.frame(x = obs_df$x, y = exp_df$y /obs_df$y )

  # I dont like this y here...

  if (monotonic){
    alt = match.arg(alternative)
    fdr$y = switch(alt,
                   left = to_monotonic(fdr$y,max=1,ascending = T),
                   right = to_monotonic(fdr$y,max=1,ascending = F)
    )
  }

  fdr$type = "fdr"

  if (any(is.infinite(fdr$y) | is.na(fdr$y))){
    fdr[ is.infinite(fdr$y) | is.na(fdr$y) ,]$y = 1 # issue #6
  }

  # this performs the substraction of the observed values against the
  # x values in the fdr dataset
  combos = abs(outer(obs,fdr$x,"-"))
  # then I get the idx corresponding to the minimum values
  min_ind = apply(combos, 1,function(x){which(x==min(x))})

  # see issue #51 for more details.
  # as I see it, this happens when a position in obs is equally similar to
  # 2 positions in the density prediction. why this only happens this
  # specific time I am unsure
  if (is.list(min_ind)){
    min_ind = unlist(lapply(min_ind, sample,size = 1))
  }

  # then I extract the corresponding fdr
  return(fdr[min_ind,]$y)
}
