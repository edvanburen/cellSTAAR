# cellSTAAR (cell-type-level STAAR)

## Citations

**cellSTAAR**:

Preprint coming soon

## Contact Emails

Eric Van Buren: [evb\@hsph.harvard.edu](mailto:evb@hsph.harvard.edu){.email}, Xihong Lin: [xlin\@hsph.harvard.edu](mailto:xlin@hsph.harvard.edu){.email}

## Introduction

<code>cellSTAAR</code> is an R package to conduct functionally informed rare variant association tests incorporating single-cell-sequencing-based functional annotations and variant sets. Given the user's own .gds file, the package can (1) create cell-type PHRED-scaled aPCs , (2) create cell-type-level variant mapping files for ENCODE cCRE categories (dELS, pELS, PLS) using each of 10 possible linking approaches, (3) run cellSTAAR and calculate cellSTAAR p-values.

## Prerequisites
<a href="https://www.r-project.org">R</a> (recommended version >= 3.5.1)

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.

## Installation

We recommend installing from GitHub or CRAN () for the latest version (1.0.2) of the package, which is built for any version of R \>= 4.0.0:

``` r
install.packages("cellSTAAR")
# OR
install.packages("devtools")
devtools::install_github("edvanburen/cellSTAAR")
```

If you are using a Mac computer and have any problems installing cellSTAAR or its required packages, one suggestion is to start by using the macrtools package (<https://github.com/coatless-mac/macrtools>) to install components that are required to compile some packages.

## cellSTAAR

cellSTAAR is summarized in the figure below: ![](/inst/image/cellSTAAR_overview.jpg)

## Usage

The workhorse function **bold** is `cellSTAAR` <code> run_cellSTAAR</code>, which can be easiest called as

-   **count_matrix**: A vector of non-negative integer counts. No normalization is done.

## Section

# subsection

## Examples

``` r
#--------------------------------------------------
#--- Simulate Data
#--------------------------------------------------

library(devtools)
devtools::install_github("edvanburen/cellSTAAR")

library(cellSTAAR)

loadRData <- function(fileName, objNameToGet = NULL){
  #loads an RData file, and returns it
  load(fileName)
  #print(ls()[ls() != "fileName"])
  if(is.null(objNameToGet)){
    rm(objNameToGet)
    #print(ls()[ls() != "fileName"])
    return(get(ls()[ls() != "fileName"]))
  }else{
    return(get(objNameToGet))
  }

}
```
## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
