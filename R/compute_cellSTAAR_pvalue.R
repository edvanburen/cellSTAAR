##' compute_cellSTAAR_pvalue.
##' @param data_obj data frame of all results from the \code{run_cellSTAAR}. By controlling the \code{grouping_vars} parameter, multiple phenotypes, cell types, linking types, and genes can be input simultaneously.
##' @param grouping_vars set of variables that uniquely identify a row in \code{data_obj}.
##' @export compute_cellSTAAR_pvalue
compute_cellSTAAR_pvalue<-function(data_obj,grouping_vars=c("gene","chr","phenotype","mapping","class","ct_name")){
  CCT_removeNA <- function(pvals, weights=NULL){
    pvals<-pvals[!is.na(pvals)]
    if(length(pvals)==0){return(NA)}
    if(!is.null(weights)){weights<-weights[!is.na(pvals)]}
    #### check if there is NA
    if(sum(is.na(pvals)) > 0){
      stop("Cannot have NAs in the p-values!")
    }

    #### check if all p-values are between 0 and 1
    if((sum(pvals<0) + sum(pvals>1)) > 0){
      stop("All p-values must be between 0 and 1!")
    }

    #### check if there are p-values that are either exactly 0 or 1.
    is.zero <- (sum(pvals==0)>=1)
    is.one <- (sum(pvals==1)>=1)
    if(is.zero && is.one){
      stop("Cannot have both 0 and 1 p-values!")
    }
    if(is.zero){
      return(0)
    }
    if(is.one){
      warning("There are p-values that are exactly 1!")
      return(1)
    }

    #### check the validity of weights (default: equal weights) and standardize them.
    if(is.null(weights)){
      weights <- rep(1/length(pvals),length(pvals))
    }else if(length(weights)!=length(pvals)){
      stop("The length of weights should be the same as that of the p-values!")
    }else if(sum(weights < 0) > 0){
      stop("All the weights must be positive!")
    }else{
      weights <- weights/sum(weights)
    }

    #### check if there are very small non-zero p-values
    is.small <- (pvals < 1e-16)
    if (sum(is.small) == 0){
      cct.stat <- sum(weights*tan((0.5-pvals)*pi))
    }else{
      cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
      cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
    }

    #### check if the test statistic is very large.
    if(cct.stat > 1e+15){
      pval <- (1/cct.stat)/pi
    }else{
      pval <- 1-pcauchy(cct.stat)
    }
    return(pval)
  }

  data_obj$pvalue<-data_obj$pval_STAAR_O
  data_obj$num_rare_var<-data_obj$num_rare_SNV


  index<-!is.na(data_obj$pval_STAAR_O_cond)
  data_obj$pvalue[index]<-data_obj$pval_STAAR_O_cond[index]
  data_obj$num_rare_var[index]<-data_obj$num_rare_SNV_cond[index]



  #t0<-data_obj%>%filter(grepl("dist",type))%>%group_by(gene,chr,phenotype,mapping,class,cutoff,ct_name)%>%mutate(CCT_pval=CCT_removeNA(pvalue),0)%>%distinct(gene,chr,phenotype,mapping,class,cutoff,ct_name,CCT_pval)%>%ungroup()
  t0<-data_obj%>%filter(grepl("dist",.data$type))%>%group_by(across(all_of(grouping_vars)))%>%mutate(CCT_pval=CCT_removeNA(.data$pvalue),0)%>%dplyr::select(all_of(grouping_vars),CCT_pval)%>%distinct()%>%ungroup()

  colnames(t0)[colnames(t0)=="CCT_pval"]<-"pvalue"
  colnames(t0)[colnames(t0)=="CCT_num_rare_var"]<-"num_rare_var"

  t0$type<-paste0("dis_combine")

  data_obj<-bind_rows(data_obj,t0)

  #t1<-data_obj%>%filter(!grepl("dist",type))%>%group_by(gene,chr,phenotype,mapping,class,cutoff,ct_name)%>%mutate(CCT_pval=CCT_removeNA(pvalue))%>%distinct(gene,chr,phenotype,mapping,class,cutoff,ct_name,CCT_pval)%>%ungroup()
  t1<-data_obj%>%filter(!grepl("dist",.data$type))%>%group_by(across(all_of(grouping_vars)))%>%mutate(CCT_pval=CCT_removeNA(.data$pvalue))%>%dplyr::select(all_of(grouping_vars),CCT_pval)%>%distinct()%>%ungroup()

  colnames(t1)[colnames(t1)=="CCT_pval"]<-"cellSTAAR_pvalue"
  t1$type<-paste0("cellSTAAR")

  return(t1)
}
