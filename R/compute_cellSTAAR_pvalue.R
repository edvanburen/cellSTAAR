##' compute_cellSTAAR_pvalue.
##' @param data_obj Data frame of all results from the \code{run_cellSTAAR} function. By controlling the \code{grouping_vars} parameter, multiple phenotypes, cell types, linking types, and genes can be input simultaneously.
##' @param grouping_vars Set of variables (Other than "link_type"!) that uniquely identify a row in \code{data_obj}.
##' @param use_conditional_p_if_exist Should the function compute the cellSTAAR omnibus p-value using conditional p-values when they exist? Defaults to FALSE.
##' @param exp_decay_dist_val Numeric value controlling the degree of exponential decay applied to the distance linking approaches for enhancer regions (pELS and dELS). If 0, all distance intervals are weighted equally. The default value of .00001 means that the six intervals 0_1, 1_50000, 50000_100000, 100000_150000, 150000_200000, and 200000_250000 are weighted with percentages of 41.4, 25.1, 15.2, 9.2, 5.6, and 3.4, respectively.
##' @export compute_cellSTAAR_pvalue

compute_cellSTAAR_pvalue<-function(data_obj,grouping_vars=c("gene","chr","phenotype","element_class","ct_name"),use_conditional_p_if_exist=FALSE,exp_decay_dist_val=.00001){
  CCT_removeNA <- function(pvals, weights=NULL){
    if(!is.null(weights)){weights<-weights[!is.na(pvals)]}
    pvals<-pvals[!is.na(pvals)]
    if(length(pvals)==0){return(NA)}
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

  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("data_obj")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }
  if("link_type"%in%grouping_vars){stop("link_type should not be included as a variable in grouping_vars input.")}

  data_obj$pvalue<-data_obj$pval_STAAR_O
  data_obj$num_rare_var<-data_obj$num_rare_SNV

  if(use_conditional_p_if_exist==TRUE){
    index<-!is.na(data_obj$pval_STAAR_O_cond)
    data_obj$pvalue[index]<-data_obj$pval_STAAR_O_cond[index]
    data_obj$num_rare_var[index]<-data_obj$num_rare_SNV_cond[index]
  }

  #browser()
  t0<-data_obj%>%filter(grepl("dist",.data$link_type))%>%
    mutate(dist_end =as.numeric(sub(".*_(\\d+)$", "\\1", link_type)))%>%
    mutate(link_weight=exp(-1*exp_decay_dist_val*dist_end))



  t0<-t0%>%group_by(across(all_of(grouping_vars)))%>%mutate(CCT_pval=CCT_removeNA(.data$pvalue,.data$link_weight),num_rare_var=round(mean(num_rare_var,na.rm=TRUE),0))%>%dplyr::select(all_of(grouping_vars),CCT_pval,num_rare_var)%>%distinct()%>%ungroup()

  colnames(t0)[colnames(t0)=="CCT_pval"]<-"pvalue"

  t0$link_type<-paste0("dis_combine")

  data_obj<-bind_rows(data_obj,t0)

  t1<-data_obj%>%filter(!grepl("dist",.data$link_type))%>%group_by(across(all_of(grouping_vars)))%>%mutate(CCT_pval=CCT_removeNA(.data$pvalue),avg_num_rare_var=round(mean(num_rare_var,na.rm=TRUE),0))%>%dplyr::select(all_of(grouping_vars),CCT_pval,avg_num_rare_var)%>%distinct()%>%ungroup()

  colnames(t1)[colnames(t1)=="CCT_pval"]<-"cellSTAAR_pvalue"
  t1$link_type<-paste0("cellSTAAR")
  return(t1)
}
