# Code to detect events.


inverse_complex_rle <- function(x, column, ...){
  if (is.null(le <- x$lengths) || is.null(v <- x[[column]]) ||
      length(le) != length(v)) {
    stop("invalid 'rle' structure")
  }
  rep.int(v, le)
}


detect_events_roberts <- function(vr,varname = "roberts_clust"){

  stopifnot(length(unique(sampleNames(vr))) == 1 )

  vr <- sortSeqlevels(vr)
  vr <- sort(vr)

  events_vector = rep(NA,length(vr))

  clust_mask = mcols(vr)[[varname]]

  events_vector[clust_mask] = "omikli"

  kat_events = clustMut::detect_events(as.numeric(!clust_mask),
                                       sig_cutoff = 0.5,
                                       event_cutoff = 5,
                                       event_value = "kataegis")

  if (length(kat_events)!=length(events_vector)){
    browser()
  }

  events_vector[!is.na(kat_events)] = kat_events[!is.na(kat_events)]

  DF = mcols(vr)
  DF$event_type = events_vector
  mcols(vr) = DF

  return(vr)
}


#' Detect events from fdr vectors
#'
#' @param x a numeric vector with lfdr values
#' @param sig_cutoff the significance cutoff
#' @param event_categories a named vector with number of the mutations in the event as values and names of the events as names.
#'
#' @return Returns a list with events and length of the stretch. Both with same total length as the fdr vector provided.
#' @export
#'
#' @examples
#'
#' fdr_test =c(runif(10,min = 0.21),
#'           c(.19,.1,.05),
#'           runif(10,min=0.21),
#'           c(.1,.1,.1,.1,.1))
#'
#' categories_test = c("kataegis" = 5,"omikli" = 2)
#'
#' detect_events(x = fdr_test,
#'   sig_cutoff = 0.2,
#'   event_categories =categories_test )
detect_events <- function(x,sig_cutoff,event_categories){


  ct_mask = x <= sig_cutoff
  event_rle = rle(ct_mask)

  event_categories = event_categories[order(event_categories)]
  events_out = rep(NA,length(event_rle$values))
  for (i in seq_along(event_categories)){
    event_name = names(event_categories)[i]
    event_cutoff = event_categories[i]
    events_out =  ifelse(event_rle$lengths >= event_cutoff & event_rle$values,
                         yes = event_name,
                         no=events_out)
  }
  ### random ids are important to then extract events afterwards
  ### random ids should be enough, we can always resplit by sampleName
  ### it's important to set the generator to FALSE to avoid giving a time
  ### dependant value
  ### I had a problem with uuids so that's why I am using radom ids which
  ### apparently uses some kind of openssl generator which don' depedent
  ### on the set.seed function. (is this ok?)
  ### length(unique(ids::random_id(1000000,bytes = 10)))
  ###  [1] 1000000
  ###  very unlikely there will be 1M of events in a sample/chr combination
  rid_vector = rep(NA,length(event_rle$values))
  rid_in = ids::random_id(n = sum(event_rle$values),
                           bytes = 10)
  rid_vector[event_rle$values] = rid_in
  event_rle$event_rid = rid_vector

  event_rle$event = events_out
  res_ev = inverse_complex_rle(event_rle,column = "event")
  res_length = inverse_complex_rle(event_rle,column = "lengths")
  res_rid = inverse_complex_rle(event_rle,column = "event_rid")
  return(list(events = res_ev,lengths = res_length,rid = res_rid))
}


