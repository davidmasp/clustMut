
<!-- README.md is generated from README.Rmd. Please edit that file -->

# clustMut <a href=''><img src='man/figures/clustmut_logo.png' align="right" height="165" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

## Overview

The goal of clustMut is to call clustered mutations in sets of somatic
mutations. A clustered mutation is defined as the result of a local
*hyper*-mutation event which generates pairs of mutations that are
closer than expected.

This **package** uses the output of a chromosomic distance randomization
from
[randommut](http://fsupeksvr.irbbarcelona.pcb.ub.es/gitlab/dmas/randommut)
and a mutation stratification list to account for subclonal frequency.

## Installation

In order to install the package use:

``` r
remotes::install_github("davidmasp/clustMut@dev")
```

You need to first install the
[remotes](https://cran.r-project.org/web/packages/remotes/index.html)
package.

You can download the bash launcher script from [this repo]() and move it
to a folder in your `$PATH`. If you update the R package, there’s no
need to update this launcher.

### Dependencies

  - Currently the package depends on genomicHelpersDMP which is beign
    transitioned to [helperMut](https://github.com/davidmasp/helperMut).
  - The rest of R packages should be installed when instaling the
    package.
  - Bash and UNIX is required to run the shell launcher script.

## Usage

A nextflow wrapper pipeline is available at
[hyperClust](https://github.com/davidmasp/hyperclust) which will contain
all the steps needed to perform the cluster calling steps.

You can use clustmut to obtain mutation clusters based on different
methods.

### distance

``` bash
clustmut distance -i /path/to/muts/ \
                --glob "*randomized.tsv" \
                --recursive \
                -o test_omichili \
                -Vlwtvu
```

### Boosting

In order to incorporate extra layers of information such as the clonal
fraction of the mutations or the strand pair a boosting list has to be
provided.

The file should be a unique column list with the following information
`{rndmut_file_name}_stratification`.

``` txt
chr:pos:sample:ref:alt_{group}
```

Where group is a set of mutations that will be grouped together in the
analysis.

NOTE: Samples should allways be included in the group as default.

## Troubleshooting

| Error code | Situation                                                                                                                                                | Solution                           |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| 123        | When the sample (or multiple samples) doesn’t have any valid mutation. This means all mutations mutations have been excluded because of the filters used | Remove that sample or the filters. |
|            |                                                                                                                                                          |                                    |

## Contribute

If you want to contribute to the package, please fork the repo and
submit a PR. Currently the package is under development so no features
are explicitily requested.
