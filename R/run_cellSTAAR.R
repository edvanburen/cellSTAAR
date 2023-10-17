##' runcellSTAAR.
##' @import STAAR
##' @import tidyverse
##' @import Matrix
##' @import gtools
##' @import doParallel
##' @import pbapply
##' @import tidyverse
##' @import dplyr
##' @import rtracklayer
##' @import SeqArray
##' @import SeqVarTools
##' @import readr
##' @importFrom stringr str_split_fixed
##' @importFrom stats quantile pcauchy pcauchy
##' @importFrom SeqArray seqOpen seqClose seqGetData
##' @importFrom SeqVarTools isSNV
##' @importFrom tibble enframe
##' @importFrom rlang .data
##' @importFrom tidyr pivot_wider
##' @importFrom utils data
##' @importFrom dplyr %>%
##' @param gds.path Path to the gds file.
##' @param ct_names ct_names.
##' @param mapping_object_list An object of class 'list' with each element being a mapping file output from the \code{create_cellSTAAR_mapping_file} function. All objects should represent the the same link approach to have structured output.
##' @param chr chromosome given as a numeric value from 1-22.
##' @param phenotype Character name of the phenotype being analyzed. Provided as part of output.
##' @param ct_aPC_list An object of class 'list' with each element being a ct_aPC object output from the \code{create_cellSTAAR_ct_aPCs} function.
##' @param null_model null model
##' @param variants_to_condition_on Dataframe of conditional variants with columns
##' @param annotation_name_catalog Dataframe with column names and locations in the .gds file for the functional annotations to include.
##' @param ncores_small Number of cores for genes with small variant sets (<500 variants)
##' @param ncores_large Number of cores for genes with large variant sets (>500 variants)
##' @param variables_to_add_to_output Dataframe of additional variables to add to output.
##' @param chr.id Used to split a chromosome into multiple jobs. Must be <= the \code{n_splits} parameter.
##' @param n_splits Number of splits for the chomosome. Used to split genes.
##' @param genes_manual Names of genes to manually analyze. Otherwise, all protein coding genes in the chromosome are analyzed.
##' @param gwas_cat_file_path gwas catalog file path.
##' @param gwas_cat_vals Values from gwas catalog file.
##' @param return_results If \code{TRUE}, the dataframe of results will be returned.
##' @param save_results If \code{TRUE}, the dataframe of results saved in the \code{out_dir} directory.
##' @param out_dir Directory to save results (used only if \code{save_results} is TRUE).
##' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01)
##' @return A data frame with the following columns, and additionally any columns passed in through the \code{variables_to_add_to_output} parameter:
##' ##' \itemize{
##' \item{\code{num_rare_SNV: }}{Number of SNV tested using \code{STAAR} function call. Equal to \code{NA} if no variants in variant set or conditional analysis run instead.}
##' \item{\code{pval_STAAR_O: }}{Omnibus p-value reported from \code{STAAR} function call. Equal to \code{NA} if no variants in variant set or conditional analysis run instead.}
##' \item{\code{num_individuals: }}{Number of individuals used in analysis.}
##' \item{\code{phenotype: }}{Phenotype analyzed.}
##' \item{\code{ct_name: }}{Name of the cell type used in variant set construction and functional annotation.}
##' \item{\code{gene: }}{Gene name.}
##' \item{\code{STAAR_time_taken: }}{Amount of time taken within \code{STAAR} function call. Equal to \code{NA} if no variants in variant set or conditional analysis run instead.}
##' ##' \item{\code{num_rare_SNV_cond: }}{Number of SNV tested using \code{STAAR_cond} function call. Equal to \code{NA} if no variants in variant set or unconditional analysis run instead.}
##' \item{\code{pval_STAAR_O_cond: }}{Omnibus p-value reported from \code{STAAR_cond} function call. Equal to \code{NA} if no variants in variant set or unconditional analysis run instead.}
##' \item{\code{n_known_var_cond: }}{Number of variants conditioned on.}
##' \item{\code{rsIDs_cond: }}{rsIDs of variants conditioned on.}
##' \item{\code{STAAR_cond_time_taken: }}{Amount of time taken within \code{STAAR_cond} function call. Equal to \code{NA} if no variants in variant set or unconditional analysis run instead.}
##' \item{\code{date: }}{Date analysis was run.}
##' \item{\code{ncores_max: }}{Maximum number of cores used in analysis.}
##' \item{\code{job_time_taken: }}{Total time taken in cellSTAAR function call.}
##' \item{\code{chr: }}{Chromosome number.}
##' }
##' @export run_cellSTAAR

run_cellSTAAR<-function(gds.path
                        ,ct_names
                       ,mapping_object_list
                       ,chr
                       ,phenotype
                        ,ct_aPC_list
                        ,null_model
                        ,variants_to_condition_on=data.frame()
                        ,annotation_name_catalog
                        ,ncores_small=1
                        ,ncores_large=1
                        ,variables_to_add_to_output=NULL
                        ,chr.id=1
                        ,n_splits=1
                        ,genes_manual=NULL
                        ,gwas_cat_file_path=NULL
                        ,gwas_cat_vals=NULL
                        ,return_results=FALSE
                        ,save_results=TRUE
                        ,out_dir="/"
                        ,rare_maf_cutoff=.01){
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("ct_names","mapping_object_list","ct_aPC_list"
                   ,"null_model","gds.path","annotation_name_catalog"
                   ,"phenotype")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if(return_results==FALSE & save_results==FALSE){
    stop("You have set both return_results and save_results as FALSE. No accessible output will be produced by the function.")
  }
  start_time<-Sys.time()
  col_names_out<-c("num_rare_SNV","pval_STAAR_O","num_individuals","phenotype","ct_name","gene","STAAR_time_taken")
  col_names_out_cond<-c("num_rare_SNV_cond","pval_STAAR_O_cond","n_known_var_cond","rsIDs_cond","ct_name","gene","STAAR_cond_time_taken")
  fit_STAAR<-function(chunk){
    #Sys.sleep(runif(1,.1,1))
    gc()
    genes_to_run<-unlist(chunk)

    results<-data.frame(matrix(NA,nrow=length(genes_to_run)*length(ct_names),ncol=length(col_names_out)))
    results_cond<-data.frame(matrix(NA,nrow=length(genes_to_run)*length(ct_names),ncol=length(col_names_out_cond)))

    results[,6]<-rep(genes_to_run,each=length(ct_names))
    results_cond[,6]<-rep(genes_to_run,each=length(ct_names))

    results[,5]<-rep(ct_names,times=length(genes_to_run))
    results_cond[,5]<-rep(ct_names,times=length(genes_to_run))

    all_pos_df2_chunk<-all_pos_df2%>%filter(.data$gene%in%genes_to_run)
    # for(g in genes_subset){
    #   all_pos_df2_chunk<-all_pos_df2%>%filter(gene%in%g)
    #   print(g)
    #   print(dim(all_pos_df2_chunk))
    # }

    chunk_unique_positions_in_use<-unique(all_pos_df2_chunk$position)
    chunk_position_index_in_use<-which(positions_in_use%in%chunk_unique_positions_in_use)
    chunk_positions_in_use<-positions_in_use[chunk_position_index_in_use]
    chunk_variantid_in_use<-variantid_in_use[chunk_position_index_in_use]
  browser()
    if(length(chunk_variantid_in_use)>0){
      for(ct_run in ct_names){
        assign(paste0("chunk_ct_aPC_",ct_run),get(paste0("ct_aPC_",ct_run))[chunk_position_index_in_use],envir = environment())
      }

      #aPC_names<-c("CADD.PHRED","LINSIGHT.PHRED","FATHMM.XF.PHRED","APC1.PHRED"
      #             ,"APC2.PHRED","APC3.PHRED","APC4.PHRED","APC5.PHRED","APC6.PHRED"
      #             ,"APC7.PHRED","APC8.PHRED","APC9.PHRED")
      aPC_names<-paste0("aPC_",anno_names)
      for(aPC_name in aPC_names)
      {
        assign(paste0("chunk_",aPC_name),get(paste0(aPC_name))[chunk_position_index_in_use])
      }
      seqSetFilter(genofile,variant.id=chunk_variantid_in_use,sample.id=pheno.id,verbose=FALSE)
      id.genotype <- seqGetData(genofile,"sample.id")
      # id.genotype.match <- rep(0,length(id.genotype))

      id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
      pheno.id.merge <- data.frame(pheno.id)
      pheno.id.merge <- dplyr::left_join(pheno.id.merge,id.genotype.merge,by=c("pheno.id"="id.genotype"))
      id.genotype.match <- pheno.id.merge$index
      ## Genotype
      chunk_Geno <- seqGetData(genofile, "$dosage")
      chunk_Geno <- chunk_Geno[id.genotype.match,]
      colnames(chunk_Geno)<-chunk_variantid_in_use

      #k<-0
      for(i in genes_to_run){
        need_conditional_analysis=FALSE
        #Min and max over any linking approach
        # since 1MB shuold not change things??
        # can tell if things changed based on what variants conditioned on anyways...
        all_pos_df2_gene<-all_pos_df2_chunk%>%filter(.data$gene==i)
        gene_unique_positions_in_use<-as.numeric(unique(all_pos_df2_gene$position))
        min_pos_set<-min(gene_unique_positions_in_use)
        max_pos_set<-max(gene_unique_positions_in_use)

        #Use forward stepwise selection procedure
        #browser()
        if(nrow(variants_to_condition_on)>0){
          cond_variant.pos<-variants_to_condition_on%>%filter(.data$POS>=min_pos_set-1e6 & .data$POS<=max_pos_set+1e6)%>%pull(.data$POS)
          cond_variant.id<-cond_var_df%>%filter(.data$cond_pos%in%cond_variant.pos)%>%pull(.data$cond_var.id)
          cond_variant.rsid<-cond_var_df%>%filter(.data$cond_pos%in%cond_variant.pos)%>%pull(.data$rsID)
        }else{
          cond_variant.pos<-numeric()
          cond_variant.id<-numeric()
          cond_variant.rsid<-numeric()
        }

        if(length(cond_variant.id)>0){
          gene_Geno_cond<-Geno_cond[,as.character(cond_variant.id),drop=FALSE]
          count_ref<-function(x){mean(x==2)}

          # Have to remove variants without any copies of the minor allele
          # to prevent an error in STAAR_cond
          var_only_ref<-apply(gene_Geno_cond,2,FUN=count_ref)
          gene_Geno_cond<-gene_Geno_cond[,var_only_ref!=1,drop=FALSE]
          n_var_adj<-NCOL(gene_Geno_cond)

          if(!is.na(n_var_adj) & n_var_adj>0){
            need_conditional_analysis=TRUE
          }
        }
        for(ct_run in ct_names){
          #k<-k+1
          all_pos_df2_gene<-all_pos_df2_chunk%>%filter(.data$gene==i,!!sym(ct_run))
          gene_unique_positions_in_use<-as.numeric(unique(all_pos_df2_gene$position))
          gene_position_index_in_use<-which(chunk_positions_in_use%in%gene_unique_positions_in_use)
          gene_positions_in_use<-chunk_positions_in_use[gene_position_index_in_use]
          gene_variantid_in_use<-chunk_variantid_in_use[gene_position_index_in_use]

          pvalues<-c()
          pvalues_cond<-c()

          difft<-NA
          cond_difft<-NA
          #n_var_adj<-NA
          #1+"e"
          if(length(gene_variantid_in_use)>=2){
            gene_ct_aPC<-get(paste0("chunk_ct_aPC_",ct_run))[gene_position_index_in_use]

            chunk_aPC_names<-paste0("chunk_",aPC_names)
            for(aPC_name in aPC_names)
            {
              assign(paste0("gene_",aPC_name),get(paste0("chunk_",aPC_name))[gene_position_index_in_use])
            }
            gene_Geno<-chunk_Geno[,as.character(gene_variantid_in_use)]


            #gene_Geno<-Geno[,as.character(gene_variantid_in_use)]

            gene_anno_cols <- lapply(paste0("gene_aPC_", anno_names), get
                                     ,envir=environment())
            anno_matrix <- data.frame(gene_ct_aPC
                                      , do.call(cbind, gene_anno_cols))
            rownames(anno_matrix)<-colnames(gene_Geno)

            min_pos_set<-min(gene_unique_positions_in_use)
            max_pos_set<-max(gene_unique_positions_in_use)

            if(!is.null(gwas_cat_file_path)){
              variants_list<-data.frame(temp_gwas2%>%dplyr::filter(.data$pos_known>=min_pos_set-1e6 & .data$pos_known<=max_pos_set+1e6)%>%arrange(.data$p_value)%>%dplyr::select(.data$CHR,.data$POS))


              # Remove any known variants in GWAS catalog from Geno
              # Should do this regardless of conditional analysis
              temp_variant.id<-gene_variantid_in_use[gene_positions_in_use%in%variants_list$POS]
              if(any(gene_variantid_in_use%in%temp_variant.id)){
                #Geno<-Geno[,-which(variant.id%in%temp_variant.id),drop=FALSE]
                gene_variantid_in_use<-gene_variantid_in_use[-which(colnames(gene_Geno)%in%temp_variant.id)]
                gene_positions_in_use<-gene_positions_in_use[-which(colnames(gene_Geno)%in%temp_variant.id)]
                gene_Geno<-gene_Geno[,-which(colnames(gene_Geno)%in%temp_variant.id),drop=FALSE]
                anno_matrix<-anno_matrix[-which(rownames(anno_matrix)%in%temp_variant.id),]
              }
            }
            #If there are any variants to potentially condition on
            # then run conditional analysis
            if(need_conditional_analysis==TRUE){
              # if TRUE, a conditional known loci is in variant set
              # and must be removed for conditional analysis
              # otherwise don't need to change Geno and anno matrices
              if(any(colnames(gene_Geno)%in%cond_variant.id)){
                gene_positions_in_use<-gene_positions_in_use[-which(colnames(gene_Geno)%in%cond_variant.id)]
                gene_variantid_in_use<-gene_variantid_in_use[-which(colnames(gene_Geno)%in%cond_variant.id)]
                Geno_nosig<-gene_Geno[,-which(colnames(gene_Geno)%in%cond_variant.id),drop=FALSE]
                anno_matrix_nosig<-anno_matrix[-which(rownames(anno_matrix)%in%cond_variant.id),]
              }else{
                Geno_nosig<-gene_Geno
                anno_matrix_nosig<-anno_matrix
              }
              print("Running Conditional STAAR")
              #1+"e"
              #stop()
              #browser()
              #print(system.time({
              st<-Sys.time()
              tryCatch({pvalues_cond <- STAAR_cond(Geno_nosig,gene_Geno_cond,null_model,anno_matrix_nosig
                                                   ,method_cond="naive",rare_maf_cutoff=rare_maf_cutoff)},error=function(e){})
              et<-Sys.time()
              cond_difft<-difftime(et,st,units="sec")
              # }))

            }
            #1+"e"
            # May not need to run this eventually
            if(!inherits(pvalues_cond,"list")){
              st<-Sys.time()
              print("Running Unconditional STAAR")
              #print(system.time({
              tryCatch({pvalues <- STAAR(gene_Geno,null_model,anno_matrix
                                         ,rare_maf_cutoff=rare_maf_cutoff)},error=function(e){})
              et<-Sys.time()
              difft<-difftime(et,st,units="sec")
              # }))
            }
          }
          # row_indexx should be equal to k, however this is a good way to ensure
          # We do not mix up the results somehow
          row_index<-which(results[,5]==ct_run&results[,6]==i)
          if(inherits(pvalues,"list")){
            results[row_index,1:4]<-c(pvalues$num_variant,pvalues$results_STAAR_O,length(pheno.id),phenotype)
            results[row_index,7]<-difft
          }else{
            results[row_index,1:4]<-c(rep(NA,2),length(pheno.id),phenotype)
            results[row_index,7]<-difft
          }
          row_index<-which(results_cond[,5]==ct_run&results_cond[,6]==i)
          if(inherits(pvalues_cond,"list")){
            results_cond[row_index,1:4]<-c(pvalues_cond$num_variant,pvalues_cond$results_STAAR_O_cond,n_var_adj,paste(cond_variant.rsid,collapse=", "))
            results_cond[row_index,7]<-cond_difft
          }else{
            results_cond[row_index,1:4]<-c(rep(NA,2),NA,NA)
            results_cond[row_index,7]<-cond_difft
          }
          #results[k,5:6]<-c(ct_run,i)
          #results_cond[k,5:6]<-c(ct_run,i)
        }

      }
    }
    gc()
    return(list(results=results,results_cond=results_cond))
  }

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

  list2env(mapping_object_list,envir = environment())

  for(ct_name in ct_names){

    gene_colnames<-colnames(get(paste0("map_obj_",ct_name)))
    genes<-gene_colnames
    genes_df<-suppressMessages(bind_cols(genes,1:length(genes)))
    colnames(genes_df)<-c("gene","col")

    index<-which(get(paste0("map_obj_",ct_name)),arr.ind=TRUE)
    index2<-suppressMessages(bind_cols(rownames(index),as_tibble(index)[-1]))
    colnames(index2)<-c("position","col")
    #browser()
    index3<-left_join(index2,genes_df,by="col")%>%dplyr::select(-.data$col)
    index3$ct<-ct_name
    assign(paste0("index_",ct_name),index3)
  }


  #print(paste("Chromosome:",chr))
  genofile <- seqOpen(gds.path)

  pheno.id <- as.character(null_model$id_include)

  id.genotype <- seqGetData(genofile,"sample.id")
  # id.genotype.match <- rep(0,length(id.genotype))

  id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
  pheno.id.merge <- data.frame(pheno.id)
  pheno.id.merge <- dplyr::left_join(pheno.id.merge,id.genotype.merge,by=c("pheno.id"="id.genotype"))%>%filter(!is.na(index))
  id.genotype.match <- pheno.id.merge$index

  ## Collect SNVs
  varid <- seqGetData(genofile, "variant.id")
  filter <- seqGetData(genofile, "annotation/filter")
  position<-seqGetData(genofile,"position")
  allele<-seqGetData(genofile,"allele")
  ALT<-gsub("*.,","",allele)
  REF<-gsub(",.*","",allele)
  isSNV<-isSNV(genofile)
  SNVlist <- filter == "PASS" & isSNV
  variant.id.SNV <- seqGetData(genofile, "variant.id")[SNVlist]
  position.SNV<-seqGetData(genofile,"position")[SNVlist]

  all_pos_df<-do.call(bind_rows,lapply(mixedsort(ls(pattern="index_",envir = environment())),get,envir=environment()))%>%distinct(.data$position,.data$gene,.data$ct)%>%mutate(value=TRUE)


  all_pos_df2<-suppressMessages(all_pos_df%>%pivot_wider(id_cols=c("position","gene")
                                                         ,names_from=.data$ct,values_from=.data$value))

  custom_fn<-function(x){sum(x,na.rm=TRUE)}

  #all_pos_df3<-all_pos_df2%>%group_by(gene)%>%dplyr::summarise(across(ct_names[1]:ct_names[length(ct_names)],\(x)sum(x,na.rm=TRUE)))

  all_pos_df3<-all_pos_df2%>%group_by(.data$gene)%>%dplyr::summarise(across(ct_names[1]:ct_names[length(ct_names)],custom_fn))


  all_pos_df3<-all_pos_df3%>%rowwise()%>%mutate(max_pos=max(across(ct_names[1]:ct_names[length(ct_names)])))%>%ungroup()

  # Put gene names in order
  # to minimize gds access & maximize computational performance

  gene_loc<-cellSTAAR::genes_biomaRt_all%>%filter(.data$gene_biotype=="protein_coding")%>%
    mutate(mid_point=.5*(.data$start_position+.data$end_position))%>%
    dplyr::select(.data$hgnc_symbol,.data$chromosome_name,.data$mid_point)%>%distinct(.data$hgnc_symbol,.keep_all = TRUE)
  colnames(gene_loc)<-c("gene","chr","mid_point")


  all_pos_df3<-left_join(all_pos_df3,gene_loc%>%dplyr::select(.data$gene,.data$mid_point),by="gene")
  gc()

  unique_positions_in_use<-unique(all_pos_df2$position)
  position_index_in_use<-which(position.SNV%in%unique_positions_in_use)
  positions_in_use<-position.SNV[position_index_in_use]
  variantid_in_use<-variant.id.SNV[position_index_in_use]
  #browser()

  if(!is.null(gwas_cat_file_path)){
    gwas_catalog<-read_tsv(gwas_cat_file_path,col_types = cols())

    temp_gwas<-gwas_catalog[gwas_catalog$`DISEASE/TRAIT`%in%gwas_cat_vals,c("DISEASE/TRAIT","CHR_ID","MAPPED_GENE","SNPS","CHR_POS","P-VALUE","DATE","CONTEXT")]
    colnames(temp_gwas)<-c("pheno","chr","mapped_gene","snp","pos_known","p_value","date","context")
    temp_gwas<-temp_gwas[!temp_gwas$chr=="X",]
    temp_gwas$p_value<-as.numeric(temp_gwas$p_value)
    temp_gwas$chr<-as.numeric(temp_gwas$chr)
    temp_gwas<-temp_gwas%>%dplyr::filter(.data$p_value<5e-8)%>%arrange(.data$chr,.data$pos_known,.data$p_value)
    #browser()
    temp_gwas<-temp_gwas[temp_gwas$chr==chr,]

    temp_gwas2<-temp_gwas
    temp_gwas2$CHR<-temp_gwas2$chr
    temp_gwas2$POS<-temp_gwas2$pos_known

    assign(paste0("temp_gwas2_",phenotype),temp_gwas2)
  }

  cond_var_df<-tibble()
  for(zzz in 1:nrow(variants_to_condition_on)){
    known_pos<-variants_to_condition_on$POS[zzz]
    known_rsid<-variants_to_condition_on$rsID[zzz]
    known_ref<-variants_to_condition_on$REF[zzz]
    known_alt<-variants_to_condition_on$ALT[zzz]

    cond_variant.id<-varid[position%in%known_pos]

    for(candidate_cond_var_id in cond_variant.id){
      r<-REF[varid%in%candidate_cond_var_id]
      a<-ALT[varid%in%candidate_cond_var_id]
      if(known_ref==r & known_alt==a){
        cond_var_df<-suppressMessages(bind_rows(cond_var_df,bind_cols(known_pos
                                                                      ,candidate_cond_var_id
                                                                      ,known_rsid)))
      }
    }
  }
  colnames(cond_var_df)<-c("cond_pos","cond_var.id","rsID")

  if(nrow(cond_var_df)>0){
    seqSetFilter(genofile,variant.id=cond_var_df$cond_var.id,sample.id=pheno.id,verbose=FALSE)
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))

    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    pheno.id.merge <- data.frame(pheno.id)
    pheno.id.merge <- dplyr::left_join(pheno.id.merge,id.genotype.merge,by=c("pheno.id"="id.genotype"))
    id.genotype.match <- pheno.id.merge$index
    Geno_cond <- seqGetData(genofile, "$dosage")
    Geno_cond <- Geno_cond[id.genotype.match,,drop=FALSE]
    colnames(Geno_cond)<-cond_var_df$cond_var.id
  }
  subgenesnum <- ceiling(length(genes)/n_splits)
  if(chr.id < n_splits)
  {
    subgenes <- ((chr.id-1)*subgenesnum + 1):(chr.id*subgenesnum)
  }
  if(chr.id == n_splits)
  {
    subgenes <- ((chr.id-1)*subgenesnum + 1):length(genes)
  }
  genes<-unique(gene_loc%>%filter(chr==chr)%>%arrange(.data$mid_point)%>%
                  filter(.data$gene%in%gene_colnames)%>%pull(.data$gene))
  genes_subset<-genes[subgenes]

  if(!is.null(genes_manual)){
    genes_subset<-genes_manual
  }
  conditional_known_loci<-vector(mode = "list", length = length(genes_subset))
  names(conditional_known_loci)<-genes_subset

  rare_maf_cutoff=.01;rv_num_cutoff=2

  j<-0
  for(ct_name in ct_names){
    j<-j+1
    #print(j)
    ct_aPC<-ct_aPC_list[[ct_name]]
    ct_aPC_in_use<-ct_aPC[position_index_in_use]
    #ct_aPC_in_use<-ct_aPC
    assign(paste0("ct_aPC_",ct_name),ct_aPC_in_use)
  }
  rm(ct_aPC_list);rm(ct_aPC);gc()
  seqSetFilter(genofile,variant.id=variantid_in_use,sample.id=pheno.id,verbose=FALSE)

  anno_names<-annotation_name_catalog$name
  #browser()
  for(i in 1:nrow(annotation_name_catalog)){
    anno_name<-annotation_name_catalog$name[i]
    anno_dir<-annotation_name_catalog$dir[i]
    assign(paste0("aPC_",anno_name),seqGetData(genofile,anno_dir))
  }

  large_genes<-all_pos_df3%>%filter(.data$max_pos>500,.data$gene%in%genes_subset)%>%arrange(.data$mid_point)%>%pull(.data$gene)
  small_genes<-setdiff(genes_subset,large_genes)

  print(paste0("Number of Small Genes: ",length(small_genes)))
  print(paste0("Number of Large Genes: ",length(large_genes)))


  #large_genes<-all_pos_df3%>%filter(max_pos>500)%>%arrange(desc(max_pos))%>%head(10)%>%pull(gene)
  small_genes<-setdiff(genes_subset,large_genes)%>%enframe()%>%dplyr::select(-.data$name)%>%dplyr::rename(gene=.data$value)
  small_genes<-left_join(small_genes,gene_loc,by=c("gene"))%>%arrange(.data$mid_point)%>%pull(.data$gene)
  #stop()
  small_size=5
  large_size=2


  small_chunks<-split(small_genes,ceiling(seq_along(small_genes)/small_size))
  large_chunks<-split(large_genes,ceiling(seq_along(large_genes)/large_size))

  n_small_chunks<-length(small_chunks)
  n_large_chunks<-length(large_chunks)
  #print(paste0("Type is ",type, " class is ", class))
  print(paste0("ncores is ",ncores_small," small genes per core is ",small_size))
  #browser()
  print(system.time({
    a<-vector('list',length=n_small_chunks)
    counter<-1
    #browser()
    while(counter<=n_small_chunks){
      out<-NULL
      start<-counter
      end<-min(counter + ncores_small - 1, n_small_chunks)
      chunks_to_run<-small_chunks[start:end]

      cl<-parallel::makeForkCluster(ncores_small)
      registerDoParallel(cl)
      #cl<-NULL
      print(paste0("Running Small Chunks ",start,"-",end," of ",n_small_chunks))
      tryCatch({out<-pblapply(chunks_to_run,FUN=fit_STAAR,cl=cl)}
               ,error=function(e){"Error in Running fit_STAAR"})
      parallel::stopCluster(cl);gc();cl=NULL
      if(is.null(out)){
        print("Chunk Failed, Trying Sequential:")
        cl<-parallel::makeForkCluster(1,outfile="out.txt")
        registerDoParallel(cl)
        tryCatch({out<-pblapply(chunks_to_run,FUN=fit_STAAR,cl=cl)}
                 ,error=function(e){print("Error in Running fit_STAAR")})
        if(is.null(out)){
          stop("Still Failed")
        }
        parallel::stopCluster(cl);gc();cl=NULL
      }
      a[start:end]<-out

      counter<-counter+ncores_small
    }
  }))
  #browser()
  #large_chunks<-large_chunks[1:8]
  print(paste0("ncores is ",ncores_large," large genes per core is ",large_size))
  print(system.time({
    b<-vector('list',length=n_large_chunks)
    counter<-1
    while(counter<=n_large_chunks){
      out<-NULL
      start<-counter
      end<-min(counter + ncores_large - 1, n_large_chunks)
      chunks_to_run<-large_chunks[start:end]

      cl<-parallel::makeForkCluster(ncores_large)
      registerDoParallel(cl)
      #cl<-NULL
      print(paste0("Running Large Chunks ",start,"-",end," of ",n_large_chunks))

      tryCatch({out<-pblapply(chunks_to_run,FUN=fit_STAAR,cl=cl)}
               ,error=function(e){"Error in Running fit_STAAR"})
      parallel::stopCluster(cl);gc();cl=NULL
      b[start:end]<-out

      counter<-counter+ncores_large
    }
  }))

  end_time<-Sys.time()

  job_time_taken<-difftime(end_time,start_time,units="secs")

  results_a<-do.call(rbind,sapply(a,'[',1))
  results_cond_a<-do.call(rbind,sapply(a,'[',2))

  results_b<-do.call(rbind,sapply(b,'[',1))
  results_cond_b<-do.call(rbind,sapply(b,'[',2))

  results<-bind_rows(results_a,results_b)
  results_cond<-bind_rows(results_cond_a,results_cond_b)

  colnames(results)<-col_names_out
  colnames(results_cond)<-col_names_out_cond
  #browser()
  results<-inner_join(results,results_cond,by=c("gene","ct_name"))
  type_noconvert_cols<-c("phenotype","gene","rsIDs_cond","ct_name")
  index<-which(!colnames(results)%in%type_noconvert_cols)

  results[,index]<-apply(results[,index],MARGIN=2,FUN=as.numeric)
  results$date<-Sys.Date()
  results$ncores_max<-ncores_small
  results$job_time_taken<-job_time_taken
  results$chr<-chr
  #browser()
  if(!is.null(variables_to_add_to_output)){
  results<-dplyr::bind_cols(results,variables_to_add_to_output)
  }
  seqClose(genofile)
  if(save_results==TRUE){
    out_name<-paste0("results_by_ct_STAAR_cell_anno")
    for(col_name in colnames(variables_to_add_to_output)){
      out_name<-paste0(out_name,"_",get(col_name))
    }
    for(ct in ct_names)
    {
      results_ct<-results%>%filter(ct_name==ct)
      out_name_ct<-paste0(out_name,"_",phenotype,"_",ct,"_chr",chr,"_",chr.id)
      assign(eval(out_name_ct),results_ct)
      save(list=eval(out_name_ct),file=paste0(out_dir,"/",out_name,".RData"))
    }
  }
  if(return_results==TRUE){
    return(results)
  }

}
