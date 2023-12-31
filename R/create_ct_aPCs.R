##' create_ct_aPCs.
##' @param gds.path Path to the gds file.
##' @param sc_epi_file_path File path to the ATAC-seq data files. It is expected that both .bw and .bed files will be in the same directory.
##' @param ct_name Name of the cell type, used for (1) loading scATAC-seq data and (2) in the file name.
##' @param num_replicate_ct_samples  Number of samples ABOVE 1. Set to NULL if the cell type has one sample, otherwise set to the total number of samples. It is expected that the samples will have similar file names: e.g. if \code{num_replicate_ct_samples=3} and \code{ct_name} is Hepatocyte, the files will have the name "Hepatocyte_1",  "Hepatocyte_2", and "Hepatocyte_3".
##' @param chr chromosome number (used as part of output filename).
##' @param out_wd Directory to save the aPCs files.
##' @return a numeric vector of cell-type PHRED-scaled values (called "aPCs" for consistency within STAAR family)
##' @export create_ct_aPCs

create_ct_aPCs<-function(gds.path
                         ,sc_epi_file_path
                         ,ct_name
                         ,num_replicate_ct_samples
                         ,chr
                         ,out_wd){
  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("gds.path","sc_epi_file_path","ct_name"
                   ,"num_replicate_ct_samples","chr","out_wd")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
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
  process_bw_aPCs<-function(path,samp_num=NULL
                            ,ct,chr_filter){
    if(!is.null(samp_num)){
      obj<-as.data.frame(import.bw(paste0(sc_epi_file_path,ct_name,"_",samp_num,".bw")))
    }else{
      obj<-as.data.frame(import.bw(paste0(sc_epi_file_path,ct_name,".bw")))
    }
    obj$seqnames<-as.character(obj$seqnames)
    obj<-obj%>%filter(.data$seqnames==chr_filter)%>%distinct(.data$seqnames,start,end,.keep_all = TRUE)
    obj$position<-0
    t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(.data$seqnames,start,end)%>%mutate(across(.data$position,~.+0:(n() - 1)))%>%ungroup()
    t1$position<-t1$position+t1$start
    t2<-t1%>%dplyr::select(.data$position,score)
    if(!is.null(samp_num)){
      colnames(t2)[2]<-paste0("score_",ct,"_",samp_num)
    }else{
      colnames(t2)[2]<-paste0("score_",ct)
    }
    t2<-right_join(t2,variant_pos_unique,by="position")
    miss<-is.na(t2[,2])
    #t2[miss,2]<-min(t2[!miss,2])/2
    t2[miss,2]<-0
    t2<-t2%>%mutate(chr=as.numeric(gsub("chr","",chr_filter)),
                    across(.cols=matches('chr'),.fns=~as.integer(.x)),
                    across(.cols=matches('position'),.fns=~as.integer(.x)))%>%dplyr::select(.data$chr,.data$position,everything())
    return(t2)
  }
  impute_anno<-function(x,val=NULL){
    if(is.null(val)){
      x[is.na(x) | x==0]<-min(x[!is.na(x) & x>0])/2
    }else{
      x[is.na(x) | x==0]<-val
    }
    return(x)
  }

  genofile <- seqOpen(gds.path)

  filter <- seqGetData(genofile, "annotation/filter")
  AVGDP <- seqGetData(genofile, "annotation/info/AVGDP")
  SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(genofile)

  variant_pos<-seqGetData(genofile,"position")[SNVlist]
  variant_pos_unique<-seqGetData(genofile,"position")[SNVlist]%>%enframe()%>%dplyr::rename(position=.data$value)%>%dplyr::select(.data$position)%>%distinct()

  # Read in ATAC-seq data
  if(!is.null(num_replicate_ct_samples)){
    j<-0
    region_file<-tibble()
    for(samp_num in 1:num_replicate_ct_samples){
      j<-j+1
      if(j==1){
        region_file<-process_bw_aPCs(path=sc_epi_file_path,samp_num=j
                                     ,ct=ct_name,chr_filter = paste0("chr",chr))
      }else{
        new<-process_bw_aPCs(path=sc_epi_file_path,samp_num=samp_num
                             ,ct=ct_name,chr_filter = paste0("chr",chr))
        region_file<-inner_join(region_file,new,by=c("chr","position"))
      }
    }

    t0<-region_file%>%mutate(score=rowMeans(across(contains(ct_name))))%>%
      dplyr::select(.data$chr,.data$position,.data$score)%>%ungroup()
    colnames(t0)[colnames(t0)=="score"]<-paste0("score_",ct_name)
  }else{
    t0<-process_bw_aPCs(path=sc_epi_file_path,ct=ct_name,chr_filter = paste0("chr",chr))
  }


    t0<-t0%>%dplyr::slice(rep(1:nrow(t0),table(variant_pos)))

    t1<-t0;t1[,3]<-impute_anno(t1[,3],val=.01)
    out1<-t1%>%pull(3);names(out1)<-NULL

    out2<-10*-log10(rank(-out1)/length(out1))

    assign(eval(paste0(ct_name,"_imputed_PHRED_chr",chr)),out2)

    out_name<-paste0(ct_name,"_imputed_PHRED_chr",chr)

    save(list=eval(out_name,envir=environment()),file=paste0(out_wd,"chr",chr,"/",out_name,".RData"))
    seqClose(genofile)
    gc()

}
