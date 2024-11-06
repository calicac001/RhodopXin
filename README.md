
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RhodopXin

## Description

`RhodopXin` is an R package designed for the exploratory analysis of
rhodopsins by looking at their sequence and structural characteristics.
Its main feature would be integrating homology modelling with motif
discovery which would help in identifying functionally important
rhodopsin-specific motifs across different species, giving an insight as
to how they have evolved across organism. Currently, there are no
pre-existing R tools that directly work with rhodopsins.

`RhodopXin` was developed using R version 4.4.1 (2024-06-14 ucrt),
Platform: x86_64-w64-mingw32/x64, and Running under: Windows 11 x64
(build 22631)

## Installation

To install the latest version of the package:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("calicac001/RhodopXin", build_vignettes = TRUE)
library("RhodopXin")
```

To run the shiny App: Under Construction

## Overview

To see an overview of the package, run:

``` r
ls("package:RhodopXin")
data(package = "RhodopXin") 
browseVignettes("RhodopXin")
```

RhodopXin contains \# functions.

1.  ***loadSequence()*** - takes in a FASTA file path and creates an
    AAStringSet object out of it

2.  ***find3dStructure()*** - given an amino acid sequence, query RCSB
    PDB to find if it has any 3D structures

3.  ***loadStructure()*** - given a rcsb id of a structure, load it in a
    viewer

4.  ***findHelices()*** - given a rcsb id of a structure, determine the
    positions of the helices

5.  ***createPWAlignments*** - generate pairwise alignments of template
    rhodopsin helices and query rhodopsins

6.  ***helixSegments*** - get a list of AAString object corresponding to
    the sequences of a rhodopsin’s helices

7.  ***findMotifs*** - (not yet implemented) identify conserved motifs
    in the rhodopsin sequence

8.  ***motifToStruct*** - (not yet implemented) map identifies motifs to
    3d structure

9.  ***basicStats*** - (not yet implemented) provide basic stats on the
    properties of the rhodopsin sequence such as length, hydrophobicity,
    and charge distribution

10. ***exploratoryPlots*** - (not yet implemented) generate exploratory
    plots such as sequence logos or motif occurrence

11. ***homology_modelling*** - (not yet implemented) allow user to
    generate 3D structural models for rhodopsin structure using a
    homology modelline if no experimental and predicted structure
    exists.

An overview of `RhodopXin` is illustrated below. TBD

## Contributions

The author of the package is Chloe Nichole Calica. The author wrote the
functions.

## References

Aboyoun P, Gentleman R (2024). *pwalign: Perform pairwise sequence
alignments*.doi:10.18129/B9.bioc.pwalign
<https://doi.org/10.18129/B9.bioc.pwalign>, R package version 1.2.0,
<https://bioconductor.org/packages/pwalign>.

Grant B.J., Rodrigues A.P.C., ElSawy K.M., McCammon J.A., Caves L.S.D.
(2006). *Bio3D: An R package for the comparative analysis of protein
structures*. <http://thegrantlab.org/bio3d/>

Korkmaz S, Yamasan B (2024). *rPDBapi: A Comprehensive Interface for
Accessing the Protein Data Bank*. R package version 2.1.1,
<https://CRAN.R-project.org/package=rPDBapi>.

Pagès H, Aboyoun P, Gentleman R, DebRoy S (2024). *Biostrings: Efficient
manipulation of biological strings*. <doi:10.18129/B9.bioc.Biostrings>
<https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.73.1,
<https://bioconductor.org/packages/Biostrings>.

Silva, A. (2022) *TestingPackage: An Example R Package For BCB410H*.
<https://github.com/anjalisilva/TestingPackage>.

Su W, Johnston B (2021). *r3dmol: Create Interactive 3D Visualizations
of Molecular Data*. R package version 0.1.2,
<https://CRAN.R-project.org/package=r3dmol>.

Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023).
*dplyr: A Grammar of Data Manipulation*. <https://dplyr.tidyverse.org>.

## Acknowledgements

This package was developed as part of an assessment for 2024 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `RhodopXin` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/calicac001/RhodopXin/issues). Many thanks to
those who provided feedback to improve this package.
