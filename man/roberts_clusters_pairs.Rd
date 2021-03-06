% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusters_alternative_methods.R
\name{roberts_clusters_pairs}
\alias{roberts_clusters_pairs}
\title{Cluster binomial for pairs}
\usage{
roberts_clusters_pairs(
  vr,
  ce_cutoff = 10,
  G_percent = 1,
  event_categories,
  dbSNP
)
}
\arguments{
\item{vr}{A VRanges object}

\item{ce_cutoff}{Complex events cutoff, only last one will be kept}

\item{G_percent}{Percentatge of the human genome used}

\item{dbSNP}{a SNPLoc object that contains the appropiate dbSNP release}
}
\value{

}
\description{
This function is an adaptations of the method in Roberts, S. A. et al. Molecular Cell (2012) to call clusterd mutations based only in pairs. Regions with local higher hypermutation than expected.
}
\details{
WARNING: This method considers the MR constant in each sample.
}
