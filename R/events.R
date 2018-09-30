# Code to detect events.


inverse_complex_rle <- function(x, column, ...){
  if (is.null(le <- x$lengths) || is.null(v <- x[[column]]) ||
      length(le) != length(v)) {
    stop("invalid 'rle' structure")
  }
  rep.int(v, le)
}

detect_events <- function(x,sig_cutoff,event_cutoff,event_value){
  ct_mask = x <= sig_cutoff
  event_rle = rle(ct_mask)
  event_rle$event = ifelse(event_rle$lengths >= event_cutoff & event_rle$values,
                           yes = event_value,
                           no=NA)
  res = inverse_complex_rle(event_rle,column = "event")
  return(res)
}

# test = runif(100)
# 
# fdr_cutoff = 0.2
# event_cutoff = 2
# type_value = "kataegis"
# detect_events(test,sig_cutoff = 0.2,event_cutoff = 2,event_value = "test")
