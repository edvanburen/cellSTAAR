# cellSTAAR (cell-type-level STAAR)

## Citations

**cellSTAAR**:

Preprint coming soon

## Contact Emails

Eric Van Buren: [evb\@hsph.harvard.edu], Xihong Lin: [xlin\@hsph.harvard.edu]

## Introduction

<code>cellSTAAR</code> is an R package to conduct functionally informed rare variant association testing incorporating single-cell-sequencing-based functional annotations and variant sets. Given an input GDS file (https://www.bioconductor.org/packages/devel/bioc/vignettes/gdsfmt/inst/doc/gdsfmt.html) and single-cell epigenetic data (in the cellSTAAR manuscript, single-cell ATAC-seq data from CATlas (http://catlas.org/humanenhancer)), the package can (1) create cell-type-level PHRED-scaled aPCs for use as functional annotation weights, (2) create cell-type-level variant mapping files for ENCODE cCRE categories (dELS, pELS, PLS) using each of 10 possible linking approaches, and (3) run cellSTAAR and calculate cellSTAAR omnibus p-values.

The current version of the <code>cellSTAAR</code> package is 1.0.1.

## Prerequisites
<a href="https://www.r-project.org">R</a> (recommended version >= 4.0.0)

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.

## Installation

We recommend installing from GitHub or CRAN () for the latest version (1.0.1) of the package, which is built for any version of R \>= 4.0.0:

``` r
#install.packages("devtools")
library(devtools)
devtools::install_github("edvanburen/cellSTAAR")
```

If you are using a Mac computer and have any problems installing cellSTAAR or its required packages, one suggestion is to start by using the macrtools package (<https://github.com/coatless-mac/macrtools>) to install components that are required to compile some packages.

Note that some dependencies for <code>cellSTAAR</code> may require installation from Bioconductor. 

## cellSTAAR

cellSTAAR is summarized in the figure below: ![](/inst/image/cellSTAAR_overview.jpg)
The key features of cellSTAAR are (1) the ability to integrate single-cell-sequencing-based functional annotations (calculated using the <code>create_ct_aPCs</code>function and variant sets (constructed using the <code>create_cellSTAAR_mapping_file</code>function) and (2) the use of the omnibus linking approach to reflect uncertainty inherent in the linking of regulatory elements to genes.
## Usage

# Create Cell-Type Variant Mapping Files
Variant mapping files for each cell type can be created using the <code>create_variant_mapping_file</code> function, which has the following input arguments:

-   **gds.path**: File path to the GDS file that will be used in the analysis
-   **sc_epi_file_path**: File path to the single-cell epigenetic files that will be used (in the manuscript, these are scATAC-seq datasets from the CATlas repository). It is expected that both .bw and .bed files will be in the same directory.
-   **ct_name**:  Name of the cell type, used for (1) loading the single-cell epigenetic data data and (2) in the created file name.
-   **num_replicate_ct_samples**: Number of samples ABOVE 1. Set to NULL if the cell type has one sample, otherwise set to the total number of samples. It is expected that the samples will have similar file names: e.g. if <code>num_replicate_ct_samples=3</code> and <code>ct_name</code> is Hepatocyte, the files will have the names "Hepatocyte_1",  "Hepatocyte_2", and "Hepatocyte_3".
-   **chr**: chromosome given as a numeric value from 1-22. This is used to filter the provided datasets and in the output name.
-   **element_class**:  One of the three ENCODE V3 cCRE categories: dELS, pELS, and PLS. Users can run dELS and pELS in one function call, but PLS must be run separately.
-  **link_types_to_run**: Character vector of one or more link types to run. The function loops over all link types specified. Comments next to the link names give the element classes for which cellSTAAR used each link. The function will throw an error if a mismatched combination of element_class and link_types_to_run is specified.
```r
c("dist_link_0_1" # pELS, dELS
   ,"dist_link_0_4000" # PLS
   ,"dist_link_1_50000" # pELS, dELS
   ,"dist_link_50000_100000" # pELS, dELS
   ,"dist_link_100000_150000" # pELS, dELS
   ,"dist_link_150000_200000" # pELS, dELS
   ,"dist_link_200000_250000" # pELS, dELS
   ,"SCREEN_link_eQTL" # pELS, dELS, PLS
   ,"SCREEN_link_noneQTL" # pELS, dELS, PLS
   ,"EpiMap_link" # pELS, dELS
   ,"ABC_link" # pELS, dELS)
````
-   **out_wd**: Directory to save the mapping files.
-   **ncores**: Number of cores to use in <code>pblapply</code> function call. Performance seems to be maximized around 3-4 cores.
-   **genes_manual**: Names of genes to manually run mapping files on. If NULL (default), all protein coding genes in the chromosome being run will be used. If specifying, ensure, the gene names used are proper HGNC symbols in the chromosome being computed.

**Cell-type-level mapping files are not phenotype specific.**

# Create Cell-Type-Level aPCs
Variant mapping files for each cell type can be created using the <code>create_ct_aPCs</code> function, which has the following input arguments:

-   **gds.path**: File path to the GDS file that will be used in the analysis
-   **sc_epi_file_path**: File path to the single-cell epigenetic files that will be used (in the manuscript, these are scATAC-seq datasets from the CATlas repository). It is expected that both .bw and .bed files will be in the same directory.
-   **ct_name**:  Name of the cell type, used for (1) loading the single-cell epigenetic data data and (2) in the created file name.
-   **num_replicate_ct_samples**: Number of samples ABOVE 1. Set to NULL if the cell type has one sample, otherwise set to the total number of samples. It is expected that the samples will have similar file names: e.g. if <code>num_replicate_ct_samples=3</code> and <code>ct_name</code> is Hepatocyte, the files will have the names "Hepatocyte_1",  "Hepatocyte_2", and "Hepatocyte_3".
-   **chr**: chromosome given as a numeric value from 1-22. This is used to filter the single-cell epigenetic datasets and in the output name.
-   **out_wd**: Directory to save the mapping files.

**Cell-type-level aPCs are not phenotype specific.**

# run cellSTAAR 
Association analysis can be run for multiple cell types simultaneously using the <code>run_cellSTAAR</code>function, which has the following input arguments:

-   **gds.path**: File path to the GDS file that will be used in the analysis
-   **ct_names**: Character vector of cell type names to run. Running multiple cell types simultaneously reduces the total computation cost by benefiting from the similarity between cell types to reduce GDS file access.
-   **chr**: Chromosome given as a numeric value from 1-22.
-   **phenotype**: Character name of the phenotype being analyzed. Provided as part of output.
-   **mapping_object_list**:  An object of class 'list' with each element being a mapping file output from the <code>create_cellSTAAR_mapping_file</code>function. All objects should represent the the same link approach to have logical output.
-   **element_class**:  One of the three ENCODE V3 cCRE categories: dELS, pELS, and PLS, corresponding to the objects input in the <code>mapping_object_list</code> argument.
-   **link_type**: Linking type name corresponding to the objects in <code>mapping_object_list </code>. See above for the acceptable values.
-   **ct_aPC_list**:  An object of class 'list' with each element being an object output from the <code>create_cellSTAAR_ct_aPCs</code>function.
-   **null_model**: Null model object output from the <code>fit_null_glmmkin</code>function of the <code>STAAR</code>package. See the examples below and the STAAR documentation (https://github.com/xihaoli/STAAR) for more details. 
-   **variants_to_condition_on**: Data frame of variants to condition on. Expected to have columns "CHR", "POS", "REF", "ALT", "rsID", and "phenotype". Defaults to an empty data frame, meaning unconditional analysis will be run for all genes. If supplied, cellSTAAR will run conditional analysis using all variants in <code>variants_to_condition_on</code>within +- 1 Mega base. 
-   **annotation_name_catalog**: Data frame with column names and locations in the GDS file for the functional annotations to include. See the examples below.
-   **analysis_to_run**: Options are "unconditional_only", "conditional_if_needed", or "both". If "conditional_if_needed", should specify <code>variants_to_condition_on</code>. Defaults to "unconditional_only" meaning no conditional analysis will be attempted. If "conditional_if_needed", only conditional analysis will be returned if there are known variants to condition on within +- 1 Mega base, otherwise only unconditional will be returned. If "both", conditional and unconditional will be attempted (note this can add substantial computation time if many genes require conditional analysis).
-  **ncores_small**: Number of cores for genes with small variant sets (<=500 variants).
-   **ncores_large**: Number of cores for genes with large variant sets (>500 variants). Larger variant sets require more memory, so most users will want to set this to be lower than <code>ncores_small </code>.
-   **variables_to_add_to_output**: Data frame of one row with additional variables to add to output. Useful for strutured output to pass into the <code>compute_cellSTAAR_pvalue</code>function.
-   **chr.id**: Used to split the genes from the analyzed chromosome into multiple jobs. Must be <= the <code>n_splits</code>parameter. Defaults to 1, meaning the entire chromosome is analyzed in one job.
-   **n_splits**: Total number of splits for genes from the chromosome being analyzed. Used to distribute computation across multiple function calls. Defaults to 1, meaning the entire chromosome is analyzed in one job.
-   **genes_manual**: Names of genes to manually run mapping files on. If NULL (default), all protein coding genes in the chromosome being run will be used. If specifying, ensure, the gene names used are proper HGNC symbols in the chromosome being computed.
-   **return_results**: If <code>TRUE</code>, the data frame of results will be returned.
-   **save_results**: If <code>TRUE</code>, the data frame of results will be saved in the <code>out_dir</code> directory.
-   **out_dir**: Directory to save results (used only if <code>save_results</code>is TRUE).
-   **rare_maf_cutoff**: The cutoff of maximum minor allele frequency in defining rare variants (default = 0.01).
-   **gwas_cat_file_path**: File path to a GWAS catalog file. This step is used to remove any rare variants that are contained within the GWAS catalog from the variant sets being tested, regardless of whether conditional analysis is used. This step is likely unnecessary, but is included to help reproduce the manuscript results. 
-   **gwas_cat_vals**: Values from the GWAS catalog corresponding to the phenotype being analyzed.

# compute cellSTAAR omnibus p-value 
The omnibus p-value from cellSTAAR can be calculated using the <code>compute_cellSTAAR_pvalue</code> function, which has the following input arguments:
-   **data_obj**: Data frame of all results from the <code>run_cellSTAAR</code>function. By controlling the <code>grouping_vars</code>parameter, multiple phenotypes, cell types, linking types, and genes can be input simultaneously.
-   **grouping_vars**: Set of variables that uniquely identify a row in <code>data_obj</code>, other than "link_type". 
-   **use_conditional_p_if_exist** Should the function compute the cellSTAAR omnibus p-value using conditional p-values when they exist? Defaults to <code>FALSE</code>.

# Examples

``` r
#-------------------------
#-- Simulate Data
#-------------------------
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
# change all chromosome values to 22
# just as an easy example
# in reality the user will likely have one GDS file
# per chromsome
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
                        
                        
#-------------------------
# Simulate cell-type-level aPCs
# Note that the function
# create_ct_aPCs is not used
# because of the small size
# of the example .gds file
#-------------------------

# Example is pretending the chromosome is 22
chr=22

# Pretending that the two cell types of interest
# are Hepatocyte and Adipocyte
ct_names<-c("Hepatocyte","Adipocyte")


ct_aPC_list<-vector('list',length=length(ct_names))
j<-0
for(ct_name in ct_names){
  j<-j+1
  ct_aPC_list[[j]]<-runif(n_variants,min=0,max=40)

}
names(ct_aPC_list)<-ct_names

#-------------------------
# Simulate variant mapping files
# Note that using the function
# create_cellSTAAR_mapping_file
# is not used because of the small size
# of the example .gds file
#-------------------------

# Take all genes from chr 22 to reflect
# file structure of a true variant mapping file
genes<-cellSTAAR::genes_biomaRt_all%>%filter(gene_biotype=="protein_coding",chromosome_name==22)%>%pull(hgnc_symbol)
n_genes<-length(genes)

# These 10 types are used for enhancers
# Promoters use dist_link_0_4000,
# SCREEN_link_eQTL, and SCREEN_link_noneQTL
link_types<-c("dist_link_0_1"
         ,"dist_link_1_50000"
         ,"dist_link_50000_100000"
         ,"dist_link_100000_150000"
         ,"dist_link_150000_200000"
         ,"dist_link_200000_250000"
         ,"SCREEN_link_eQTL"
         ,"SCREEN_link_noneQTL"
         ,"EpiMap_link"
         ,"ABC_link")
         
n_cts<-length(ct_names)

for(link_type in link_types){
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
  assign(paste0("map_objs_",link_type),map_objs)
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

#-------------------------
# Run all linking approaches
#-------------------------

# Note that these names correspond to the
# mapping object lists created above

link_types<-c("dist_link_0_1"
        ,"dist_link_1_50000"
        ,"dist_link_50000_100000"
        ,"dist_link_100000_150000"
        ,"dist_link_150000_200000"
        ,"dist_link_200000_250000"
        ,"SCREEN_link_eQTL"
        ,"SCREEN_link_noneQTL"
        ,"EpiMap_link"
        ,"ABC_link")
         
#Counter of across types loop    
j<-0
for(link_type in link_types){
  j<-j+1
  # Will be added to output to aid in computing
  # cellSTAAR omnibus p-value using the
  # compute_cellSTAAR_pvalue function
  variable_df<-dplyr::bind_rows(element_source="cCRE_V3"
                                ,"sc_cutoff"="0.8")
  
  print(paste0("Type ",link_type,"; # ", j, " of ",length(link_types)))
  assign(paste0("results_cellSTAAR_",link_type),run_cellSTAAR(ct_names
                                         ,genes_manual=genes[1:5] #run five genes as an example
                                         ,chr=chr
                                         ,phenotype = "PHENO"
                                         ,mapping_object_list=get(paste0("map_objs_",link_type))
                                         ,link_type=link_type
                                         ,element_class="dELS"
                                         ,ct_aPC_list=ct_aPC_list
                                         ,null_model=null_model
                                         ,gds.path=gds.path
                                         ,variants_to_condition_on=data.frame() # Unconditional analysis
                                         ,annotation_name_catalog=annotation_name_catalog
                                         ,variables_to_add_to_output=variable_df
                                         ,save_results = FALSE #do not save to disk
                                         ,return_results = TRUE #rather return object since quick example
                                         ,rare_maf_cutoff=1 # Just a syntax example, cellSTAAR should be 
                                                            # used to test only rare variants
                                         ))
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
