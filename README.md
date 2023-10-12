---
output:
  html_document: default
  pdf_document: default
---
# cellSTAAR

## Citations

**cellSTAAR**:

Preprint coming soon

## Contact Emails

Eric Van Buren: [evb\@hsph.harvard.edu](mailto:evb@hsph.harvard.edu){.email}, Xihong Lin: [xlin\@hsph.harvard.edu](mailto:xlin@hsph.harvard.edu){.email}

## Introduction

<code>cellSTAAR</code> is an R package to conduct functionally informed rare variant association tests incorporating single-cell-sequencing-based functional annotations and variant sets.

## Installation

We recommend installing from GitHub or CRAN () for the latest version (1.0.2) of the package, which is built for any version of R \>= 4.0.0:

``` r
install.packages("cellSTAAR")
# OR
install.packages("devtools")
devtools::install_github("edvanburen/cellSTAAR")
```

Note the following minimum package versions for imported packages: multcomp (\>= 1.4-13), methods, pscl (\>= 1.5.5), pbapply (\>= 1.4.0), parallel (\>= 3.6.3), doParallel (\>= 1.0.15).

If you are using a Mac computer and have any problems installing cellSTAAR or its required packages, one suggestion is to start by using the macrtools package (<https://github.com/coatless-mac/macrtools>) to install components that are required to compile some packages.

## cellSTAAR

cellSTAAR is summarized in the figure below: ![](/inst/image/cellSTAAR_overview.jpg){align="center"}

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
```
