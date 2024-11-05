
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RhodopXin

## Description

`RhodopXin` is an R package designed for the exploratory analysis of
rhodopsins by looking at their structural, functional, and sequence
characteristics. Its main feature would be integrating homology
modelling with motif discovery which would help in identifying
functionally important rhodopsin-specific motifs across different
species, giving an insight as to how they have evolved across organism.

Currently, there are no pre-existing R tools that directly work with
rhodopsins. Moreover, this package would leverage pre-existing
bioinformatics tools for protein analysis and extend their usage to
something that is more specific i.e. exploring the properties of
rhodopsins in a comprehensive manner.

## Installation

To install the latest version of the package:

``` r
install.packages("devtools")
#> Installing package into 'C:/Users/chloe/AppData/Local/Temp/RtmpERKag8/temp_libpath338576e2045'
#> (as 'lib' is unspecified)
#> package 'devtools' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\chloe\AppData\Local\Temp\RtmpWU43vV\downloaded_packages
library("devtools")
#> Loading required package: usethis
devtools::install_github("calicac001/RhodopXin", build_vignettes = TRUE)
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo calicac001/RhodopXin@HEAD
#> pkgbuild  (1.4.4  -> 1.4.5   ) [CRAN]
#> fs        (1.6.4  -> 1.6.5   ) [CRAN]
#> glue      (1.7.0  -> 1.8.0   ) [CRAN]
#> withr     (3.0.1  -> 3.0.2   ) [CRAN]
#> waldo     (0.5.3  -> 0.6.0   ) [CRAN]
#> rlang     (1.1.3  -> 1.1.4   ) [CRAN]
#> ps        (1.8.0  -> 1.8.1   ) [CRAN]
#> jsonlite  (1.8.8  -> 1.8.9   ) [CRAN]
#> evaluate  (1.0.0  -> 1.0.1   ) [CRAN]
#> cli       (3.6.2  -> 3.6.3   ) [CRAN]
#> Rcpp      (1.0.13 -> 1.0.13-1) [CRAN]
#> sys       (3.4.2  -> 3.4.3   ) [CRAN]
#> askpass   (1.2.0  -> 1.2.1   ) [CRAN]
#> openssl   (2.2.1  -> 2.2.2   ) [CRAN]
#> curl      (5.2.2  -> 6.0.0   ) [CRAN]
#> tinytex   (0.52   -> 0.54    ) [CRAN]
#> xfun      (0.45   -> 0.49    ) [CRAN]
#> rmarkdown (2.28   -> 2.29    ) [CRAN]
#> bio3d     (2.4-4  -> 2.4-5   ) [CRAN]
#> Skipping 2 packages not available: msa, Biostrings
#> Installing 19 packages: pkgbuild, fs, glue, withr, waldo, rlang, ps, jsonlite, evaluate, cli, Rcpp, sys, askpass, openssl, curl, tinytex, xfun, rmarkdown, bio3d
#> Installing packages into 'C:/Users/chloe/AppData/Local/Temp/RtmpERKag8/temp_libpath338576e2045'
#> (as 'lib' is unspecified)
#> 
#>   There are binary versions available but the source versions are later:
#>           binary source needs_compilation
#> waldo      0.5.3  0.6.0             FALSE
#> curl       5.2.3  6.0.0              TRUE
#> rmarkdown   2.28   2.29             FALSE
#> 
#> package 'pkgbuild' successfully unpacked and MD5 sums checked
#> package 'fs' successfully unpacked and MD5 sums checked
#> package 'glue' successfully unpacked and MD5 sums checked
#> package 'withr' successfully unpacked and MD5 sums checked
#> package 'rlang' successfully unpacked and MD5 sums checked
#> package 'ps' successfully unpacked and MD5 sums checked
#> package 'jsonlite' successfully unpacked and MD5 sums checked
#> package 'evaluate' successfully unpacked and MD5 sums checked
#> package 'cli' successfully unpacked and MD5 sums checked
#> package 'Rcpp' successfully unpacked and MD5 sums checked
#> package 'sys' successfully unpacked and MD5 sums checked
#> package 'askpass' successfully unpacked and MD5 sums checked
#> package 'openssl' successfully unpacked and MD5 sums checked
#> package 'tinytex' successfully unpacked and MD5 sums checked
#> package 'xfun' successfully unpacked and MD5 sums checked
#> package 'bio3d' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\chloe\AppData\Local\Temp\RtmpWU43vV\downloaded_packages
#> installing the source packages 'waldo', 'curl', 'rmarkdown'
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\chloe\AppData\Local\Temp\RtmpWU43vV\remotes43ec7d3a5aed\calicac001-RhodopXin-333552c/DESCRIPTION' ...  ✔  checking for file 'C:\Users\chloe\AppData\Local\Temp\RtmpWU43vV\remotes43ec7d3a5aed\calicac001-RhodopXin-333552c/DESCRIPTION' (341ms)
#>       ─  preparing 'RhodopXin':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>      Omitted 'LazyData' from DESCRIPTION
#>       ─  building 'RhodopXin_0.1.0.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/chloe/AppData/Local/Temp/RtmpERKag8/temp_libpath338576e2045'
#> (as 'lib' is unspecified)
library("RhodopXin")
```

To run the shiny App: Under Construction

## Overview

To see an overview of the package, run:

``` r
ls("package:RhodopXin")
#> [1] "loadSequence"
data(package = "RhodopXin") 
#> no data sets found
browseVignettes("RhodopXin")
#> No vignettes found by browseVignettes("RhodopXin")
```

RhodopXin contains \# functions.

1.  

An overview of `RhodopXin` is illustrated below. TBD

## Contributions

The author of the package is Chloe Nichole Calica. The author wrote the
functions.

## References

Silva, A. (2022) TestingPackage: An Example R Package For BCB410H.
<https://github.com/anjalisilva/TestingPackage>

Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023).
dplyr: A Grammar of Data Manipulation. <https://dplyr.tidyverse.org>

## Acknowledgements

This package was developed as part of an assessment for 2024 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `RhodopXin` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/calicac001/RhodopXin/issues). Many thanks to
those who provided feedback to improve this package.
