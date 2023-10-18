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

We recommend installing from GitHub or CRAN () for the latest version (1.0.1) of the package, which is built for any version of R \>= 4.0.0:

``` r
#install.packages("devtools")
library(devtools)
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
# Load required packages to simulate example
# .gds data file, kinship matrix, null model,
# and cell-type aPCs
library(cellSTAAR)
library(SeqVarTools)
library(dplyr)
library(Matrix)
library(kinship2)

#Set seed for reproducibility
set.seed(1234)

# Use example data in SeqVarTools package
# to create a sample .gds file
# and open it
vcffile <- seqExampleFileName("vcf")
gds.path <- "tmp.gds"
seqVCF2GDS(vcffile, gds.path, verbose=FALSE)
gds<-seqOpen(gds.path,readonly = FALSE)

# Collect information from sample .gds file
pos<-seqGetData(gds,"position")
n_variants<-length(pos)
# 90 individuals
sample.id<-seqGetData(gds,"sample.id")
variant.id<-seqGetData(gds,"variant.id")
n_ind<-length(sample.id)

# Simulate average depth and 3 aPCs 
AVGDP<-rep(20,n_variants)
aPC1<-runif(n_variants,min=0,max=40)
aPC2<-runif(n_variants,min=0,max=40)
aPC3<-runif(n_variants,min=0,max=40)

# Add simulated data into .gds file
# while pretending the chromosome is 22
Anno.folder <- index.gdsn(gds, "annotation/info")
add.gdsn(Anno.folder, "AVGDP", val=AVGDP, compress="")
add.gdsn(Anno.folder, "aPC1", val=aPC1, compress="")
add.gdsn(Anno.folder, "chr", val=rep(22,n_variants), compress="")
add.gdsn(Anno.folder, "aPC2", val=aPC2, compress="")
add.gdsn(Anno.folder, "aPC3", val=aPC3, compress="")
seqClose(gds)

# Simulate a phenotype ("PHENO") 
# and three covariates for null model
PHENO<-rnorm(n_ind,mean=10,sd=1)
PC1<-rnorm(n_ind,mean=0,sd=1)
PC2<-rnorm(n_ind,mean=0,sd=1)
sex<-rbinom(n_ind,1,.5)
pheno_data<-bind_cols(sample.id,PHENO,PC1,PC2,sex)
colnames(pheno_data)<-c("sample.id","PHENO","PC1","PC2","sex")

# Create a kinship matrix representing
# 15 families of 6 individuals each
# this choice was completely arbitrary
# and designed to correspond to the 90
# individuals in the sample VCF/GDS file
grid <- 1
Npercell <- n_ind
ndiv <- 1
vfam <- 0.5
N <- round(grid*grid*Npercell/ndiv)

unitmat <- matrix(0.5, 4, 4)
diag(unitmat) <- 1
unitmat[1,2] <- unitmat[2,1] <- 0
ped <- data.frame(famid = rep(as.integer(1:15), each=6), id = as.integer(1:n_ind), fa = rep(0, n_ind), mo = rep(0, n_ind))
for(i in 1:15) {
  ped$fa[6*i-(0:1)] <- ped$id[6*i-3]
  ped$mo[6*i-(0:1)] <- ped$id[6*i-2]
}
kins <- makekinship(ped$famid, ped$id, ped$fa, ped$mo)
rownames(kins)<-sample.id
colnames(kins)<-sample.id

#Fit the null model using the STAAR package
null_model<-STAAR::fit_null_glmmkin(PHENO~PC1+PC2+sex,use_sparse=TRUE
                        ,kins=kins,data=pheno_data,id="sample.id"
                        ,family=gaussian(link="identity"))
                        
                        
####################################
# Simulate cell-type-level aPCs
# Note that the function
# create_ct_aPCs is not used
# because of the small size
# of the example .gds file
####################################
chr=22
ct_names<-c("Hepatocyte","Adipocyte")
ct_aPC_list<-vector('list',length=length(ct_names))
j<-0
for(ct_name in ct_names){
  j<-j+1
  ct_aPC_list[[j]]<-runif(n_variants,min=0,max=40)

}
names(ct_aPC_list)<-ct_names

####################################
# Simulate variant mapping files
# Note that using the function
# create_cellSTAAR_mapping_file
# is not used because of the small size
# of the example .gds file
####################################

# Take all genes from chr 22 to better reflect
# file structure of a true variant mapping file
genes<-cellSTAAR::genes_biomaRt_all%>%filter(gene_biotype=="protein_coding",chromosome_name==22)%>%pull(hgnc_symbol)
n_genes<-length(genes)

types<-c("dist_0_1_filter_CATlas","dist_1_50000_filter_CATlas"
         ,"dist_50000_100000_filter_CATlas","dist_100000_150000_filter_CATlas"
         ,"dist_150000_200000_filter_CATlas"
         ,"dist_200000_250000_filter_CATlas","ABC_link_filter_CATlas"
         ,"EpiMap_link_filter_CATlas","SCREEN_link_eQTL_filter_CATlas"
         ,"SCREEN_link_noneQTL_filter_CATlas")
         
n_cts<-length(ct_names)

for(type in types){
  map_objs<-vector('list',length=n_cts)
  j<-0
  for(ct_name in ct_names){
    j<-j+1
    mat<-Matrix(matrix(as.logical(rbinom(n_variants * n_genes, 1, 0.1)), n_variants, n_genes),sparse=TRUE)
    colnames(mat)<-genes
    rownames(mat)<-pos
    map_objs[[j]]<-mat
    names(map_objs)[j]<-paste0("map_obj_",ct_name)
  }
  assign(paste0("map_objs_",type),map_objs)
  rm(map_objs)
}

# Set up "annotation_name_catalog"
# This gives the names and locations
# of annotations inside the .gds file
# to use as weights in the STAAR procedure
annotation_names<-c("aPC1","aPC2","aPC3")
annotation_locations<-c("annotation/info/aPC1","annotation/info/aPC2","annotation/info/aPC3")
annotation_name_catalog<-dplyr::bind_cols(annotation_names,annotation_locations)
colnames(annotation_name_catalog)<-c("name","dir")

# Run all types
types<-c("dist_0_1_filter_CATlas","dist_1_50000_filter_CATlas"
         ,"dist_50000_100000_filter_CATlas","dist_100000_150000_filter_CATlas"
         ,"dist_150000_200000_filter_CATlas"
         ,"dist_200000_250000_filter_CATlas","ABC_link_filter_CATlas"
         ,"EpiMap_link_filter_CATlas","SCREEN_link_eQTL_filter_CATlas"
         ,"SCREEN_link_noneQTL_filter_CATlas")



gwas_cat_file_path="/n/holystore01/LABS/xlin/Lab/evb/data/gwas_catalog_v1.0-associations_e100_r2021-02-25.tsv"
j<-0
for(type in types){
  variable_df<-dplyr::bind_rows("mapping"="cCRE_V3"
                                ,"class"="pretend_dELS" 
                                ,"type"=type
                                ,"cutoff"="0.8")
  j<-j+1
  print(paste0("Type ",type,"; # ", j, " of ",length(types)))
  assign(paste0("results_cellSTAAR_",type),run_cellSTAAR(ct_names
                                         ,genes_manual=genes[1:5] #run five genes as an example
                                         ,chr=chr
                                         ,phenotype = "PHENO"
                                         ,mapping_object_list=get(paste0("map_objs_",type))
                                         ,ct_aPC_list=ct_aPC_list
                                         ,null_model=null_model
                                         ,gds.path=gds.path
                                         ,variants_to_condition_on=data.frame() # Unconditional analysis
                                         ,annotation_name_catalog=annotation_name_catalog
                                         ,variables_to_add_to_output=variable_df
                                         ,save_results = FALSE #do not save to disk
                                         ,return_results = TRUE #rather return object since quick example
                                         ,rare_maf_cutoff=1))
  #rm(list=paste0("map_objs_",type))
  gc()
}

# Add all results into one data frame
# and compute cellSTAAR p-value
assign(paste0("all"),do.call(dplyr::bind_rows,lapply(gtools::mixedsort(ls(pattern="results_",envir = environment())),get,envir=environment())))

cellSTAAR_p_values<-compute_cellSTAAR_pvalue(all)

all<-bind_rows(all,cellSTAAR_p_values)

```
## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
