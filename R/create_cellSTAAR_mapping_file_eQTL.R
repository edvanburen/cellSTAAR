##' create_cellSTAAR_mapping_file.
##' @param gds.path Path to the gds file.
##' @param eQTL_match_file match file for eQTL summary statistics.
##' @param chr chromosome given as a numeric value from 1-22.
##' @param element_class One of the three ENCODE V3 cCRE categories: dELS, pELS, and PLS. Users can run dELS and pELS in one function call, but PLS must be run separately.
##' @param link_types_to_run Character vector of link types to run. The function loops over all link types specified.
##' @param out_wd Directory to save the mapping files.
##' @param ncores Number of cores to use in \code{pblapply} call.
##' @param genes_manual Names of genes to manually run mapping files on. If NULL, all protein coding genes in the chromosome being run will be used. If specifying, ensure, the gene names used are proper HGNC symbols in the chromosome being computed.

##' @return a sparse matrix, with rows covering variant positions and colnames covering protein coding genes. A value of 1 indicates a link between the respective position and the particular gene, 0 indicates no link.
##' @export create_cellSTAAR_mapping_file_eQTL
create_cellSTAAR_mapping_file_eQTL<-function(gds.path
                                             ,eQTL_match_file
                                        ,chr
                                        ,element_class
                                        ,link_types_to_run
                                        ,out_wd
                                        ,ncores=1
                                        ,genes_manual=NULL){

  passed_args <- names(as.list(match.call())[-1])
  required_args<-c("gds.path"
                   ,"chr","link_types_to_run","element_class"
                   ,"out_wd")
  if (any(!required_args %in% passed_args)) {
    stop(paste("Argument(s)",paste(setdiff(required_args, passed_args), collapse=", "),"missing and must be specified."))
  }

  for(ec in element_class){
    if(!ec%in%c("dELS","pELS","PLS")){
      stop(paste0("element class must be either dELS, pELS, or PLS"))}
  }


  for(lt_to_check in link_types_to_run){
    if(!lt_to_check%in%c("dist_link_0_1"
                         ,"dist_link_0_4000"
                         ,"dist_link_1_50000"
                         ,"dist_link_50000_100000"
                         ,"dist_link_100000_150000"
                         ,"dist_link_150000_200000"
                         ,"dist_link_200000_250000"
                         ,"SCREEN_link_eQTL"
                         ,"SCREEN_link_noneQTL"
                         ,"EpiMap_link"
                         ,"ABC_link")){
      stop(paste0("link_type must be one of ",paste0(c("dist_link_0_1"
                                                       ,"dist_link_0_4000"
                                                       ,"dist_link_1_50000"
                                                       ,"dist_link_50000_100000"
                                                       ,"dist_link_100000_150000"
                                                       ,"dist_link_150000_200000"
                                                       ,"dist_link_200000_250000"
                                                       ,"SCREEN_link_eQTL"
                                                       ,"SCREEN_link_noneQTL"
                                                       ,"EpiMap_link"
                                                       ,"ABC_link"),collapse=" ")))
    }
    if(grepl("pELS|dELS",element_class) & lt_to_check%in%c("dist_link_0_4000")){
      stop("Link type should not be dist_link_0_4000 when consructing mapping files for element class pELS or dELS. Please call the function separately for enhancers and promoters.")
    }
    if(grepl("PLS",element_class) & lt_to_check%in%c("dist_link_0_1"
                                                    ,"dist_link_1_50000"
                                                    ,"dist_link_50000_100000"
                                                    ,"dist_link_100000_150000"
                                                    ,"dist_link_150000_200000"
                                                    ,"dist_link_200000_250000"
                                                    ,"EpiMap_link"
                                                    ,"ABC_link")){
      stop("Link type should only be dist_link_0_4000, SCREEN_link_eQTL, or SCREEN_link_non eQTL when consructing mapping files for element class PLS. Please call the function separately for enhancers and promoters.")
    }
  }


  genofile <- seqOpen(gds.path)

  filter <- seqGetData(genofile, "annotation/filter")
  AVGDP <- seqGetData(genofile, "annotation/info/AVGDP")
  SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(genofile)
  positions<-seqGetData(genofile,"position")[SNVlist]
  variant_pos<-positions%>%enframe()%>%dplyr::rename(position=.data$value)%>%dplyr::select(.data$position)%>%distinct()

  for(link_type in link_types_to_run){
    if(link_type=="SCREEN_link_eQTL"){
      raw_mappings_SCREEN<-cellSTAAR::agnostic_dnase_summary_V3_eQTL%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene,.keep_all = TRUE)%>%filter(.data$gene!="")
    }
    if(link_type=="SCREEN_link_noneQTL"){
      raw_mappings_SCREEN<-cellSTAAR::agnostic_dnase_summary_V3_noneQTL%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene,.keep_all = TRUE)%>%filter(.data$gene!="")
    }
    if(grepl("dist_link_0_4000",link_type)){
      raw_mappings_dist<-cellSTAAR::raw_mappings_cCRE_V3_dist_0_4000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene_dist_0_4000,.keep_all = TRUE)
    }
    if(grepl("dist_link_0_1",link_type)){
      raw_mappings_dist<-cellSTAAR::raw_mappings_cCRE_V3_dist_0_1%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene_dist_0_1,.keep_all = TRUE)
    }
    if(grepl("dist_link_1_50000",link_type)){
      raw_mappings_dist<-cellSTAAR::raw_mappings_cCRE_V3_dist_1_50000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene_dist_1_50000,.keep_all = TRUE)
    }
    if(grepl("dist_link_50000_100000",link_type)){
      raw_mappings_dist<-cellSTAAR::raw_mappings_cCRE_V3_dist_50000_100000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene_dist_50000_100000,.keep_all = TRUE)
    }
    if(grepl("dist_link_100000_150000",link_type)){
      raw_mappings_dist<-cellSTAAR::raw_mappings_cCRE_V3_dist_100000_150000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene_dist_100000_150000,.keep_all = TRUE)
    }
    if(grepl("dist_link_150000_200000",link_type)){
      raw_mappings_dist<-cellSTAAR::raw_mappings_cCRE_V3_dist_150000_200000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene_dist_150000_200000,.keep_all = TRUE)
    }
    if(grepl("dist_link_200000_250000",link_type)){
      raw_mappings_dist<-cellSTAAR::raw_mappings_cCRE_V3_dist_200000_250000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,.data$cCRE_accession,.data$gene_dist_200000_250000,.keep_all = TRUE)
    }
    if(link_type=="EpiMap_link"){
      raw_mappings_EpiMap<-cellSTAAR::raw_mappings_cCRE_V3_EpiMap_link_all_50%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,.data$cCRE_accession,.data$EpiMap_gene,.keep_all = TRUE)
    }
    if(link_type=="ABC_link"){
      raw_mappings_ABC<-cellSTAAR::raw_mappings_cCRE_V3_ABC_link_all_50%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,.data$cCRE_accession,.data$ABC_gene,.keep_all = TRUE)
    }

    split_into_multiple <- function(column, pattern = ", ", into_prefix){
      cols <- str_split_fixed(column, pattern, n = Inf)
      # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
      cols[which(cols == "")] <- NA
      cols <- as_tibble(cols)
      # name the 'cols' tibble as 'into_prefix_1', 'into_prefix_2', ..., 'into_prefix_m'
      # where m = # columns of 'cols'
      m <- dim(cols)[2]

      names(cols) <- paste(into_prefix, 1:m, sep = "")
      return(cols)
    }
    bp_level_mappings<-function(obj,filt,other_cts=FALSE){
      #This is protection in case a group tibble is passed in
      # Running this code on a grouped tibble will cause silent errors!!
      obj<-ungroup(obj)
      if(nrow(obj)==0){return(obj)}
      obj$position<-0

      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(.data$cCRE_accession,.data$gene)%>%mutate(across(.data$position,~.+0:(n() - 1)))%>%ungroup()
      t1$position<-t1$position+t1$start

      t2<-t1%>%arrange(.data$position)%>%dplyr::select(.data$position,.data$gene,.data$classification1,.data$cCRE_accession)%>%group_by(.data$position)%>%mutate(genes=paste(.data$gene,collapse=","))%>%dplyr::select(.data$position,.data$genes,.data$classification1,.data$cCRE_accession)%>%arrange(.data$position)%>%distinct(.data$position,.data$genes,.data$classification1,.data$cCRE_accession)
      t3<-t2%>%bind_cols(split_into_multiple(column=.$genes,pattern=",",into_prefix="cCRE_gene"))
      colnames(t3)[2]<-paste0("genes_cCRE")

      t3<-inner_join(t3,variant_pos,by="position")%>%distinct()
      if(filt=="CATlas"){
        # Keeps any CATlas peak region
        t5<-inner_join(ct_CATlas_pos,t3,by="position")%>%distinct()

      }
      if(filt%in%c("nofilter")){
        t5<-t3
      }
      return(t5)
    }
    bp_level_mappings_dist<-function(obj,filt,other_cts=FALSE){
      #This is protection in case a group tibble is passed in
      # Running this code on a grouped tibble will cause silent errors!!
      obj<-ungroup(obj)
      if(nrow(obj)==0){return(obj)}
      obj$position<-0
      dist_val<-gsub("dist_link_","",link_type)
      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(.data$cCRE_accession,!!as.symbol(paste0("gene_dist_",dist_val)))%>%mutate(across(.data$position,~.+0:(n() - 1)))%>%ungroup()
      t1$position<-t1$position+t1$start

      t2<-t1%>%arrange(.data$position)%>%dplyr::select(.data$position,!!as.symbol(paste0("gene_dist_",dist_val)),.data$classification1,.data$cCRE_accession)%>%distinct()
      colnames(t2)[2]<-"cCRE_gene1"
      #t3<-t2%>%bind_cols(split_into_multiple(column=.$genes,pattern=",",into_prefix="cCRE_gene"))
      #colnames(t3)[2]<-paste0("genes_cCRE")
      t3<-t2
      t3<-inner_join(t3,variant_pos,by="position")%>%distinct()

      if(filt=="CATlas"){
        # Keeps any CATlas peak region
        t5<-inner_join(ct_CATlas_pos,t3,by="position")%>%distinct()

      }
      if(filt%in%c("nofilter")){
        t5<-t3
      }
      return(t5)
    }
    bp_level_mappings_EpiMap_link<-function(obj,filt,other_cts=FALSE){
      #This is protection in case a group tibble is passed in
      # Running this code on a grouped tibble will cause silent errors!!
      obj<-ungroup(obj)
      if(nrow(obj)==0){return(obj)}
      obj$position<-0

      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(.data$cCRE_accession,.data$EpiMap_gene)%>%mutate(across(.data$position,~.+0:(n() - 1)))%>%ungroup()

      t1$position<-t1$position+t1$start

      ##If a position is in multiple regions we really need to remove one
      #temp_to_keep<-t1%>%arrange(.data$position,.data$cCRE_accession)%>%distinct(.data$position,.data$cCRE_accession)%>%group_by(.data$position)%>%dplyr::slice(1)
      #t1<-right_join(t1,temp_to_keep,by=c(".data$cCRE_accession","position"))
      t2<-t1%>%arrange(.data$position)%>%dplyr::select(.data$position,.data$EpiMap_gene,.data$classification1,.data$cCRE_accession)%>%group_by(.data$position)%>%mutate(genes=paste(.data$EpiMap_gene,collapse=","))%>%dplyr::select(.data$position,.data$genes,.data$classification1,.data$cCRE_accession)%>%arrange(.data$position)%>%distinct(.data$position,.data$genes,.data$classification1,.data$cCRE_accession)
      t3<-t2%>%bind_cols(split_into_multiple(column=.$genes,pattern=",",into_prefix="EpiMap_gene"))
      colnames(t3)[2]<-paste0("genes_EpiMap")

      #t3<-t2
      t3<-inner_join(t3,variant_pos,by="position")%>%distinct()


      if(filt=="CATlas"){
        # Keeps any CATlas peak region
        t5<-inner_join(ct_CATlas_pos,t3,by="position")%>%distinct()

      }
      if(filt%in%c("nofilter")){
        t5<-t3
      }
      return(t5)
    }
    bp_level_mappings_ABC_link<-function(obj,filt,other_cts=FALSE){
      #This is protection in case a group tibble is passed in
      # Running this code on a grouped tibble will cause silent errors!!
      obj<-ungroup(obj)
      if(nrow(obj)==0){return(obj)}
      obj$position<-0

      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(.data$cCRE_accession,.data$ABC_gene)%>%mutate(across(.data$position,~.+0:(n() - 1)))%>%ungroup()
      t1$position<-t1$position+t1$start

      t2<-t1%>%arrange(.data$position)%>%dplyr::select(.data$position,.data$ABC_gene,.data$classification1,.data$cCRE_accession)%>%group_by(.data$position)%>%mutate(genes=paste(.data$ABC_gene,collapse=","))%>%dplyr::select(.data$position,.data$genes,.data$classification1,.data$cCRE_accession)%>%arrange(.data$position)%>%distinct(.data$position,.data$genes,.data$classification1,.data$cCRE_accession)
      t3<-t2%>%bind_cols(split_into_multiple(column=.$genes,pattern=",",into_prefix="ABC_gene"))

      colnames(t3)[2]<-paste0("genes_ABC")
      t3<-inner_join(t3,variant_pos,by="position")%>%distinct()


      if(filt=="CATlas"){
        # Keeps any CATlas peak region
        t5<-inner_join(ct_CATlas_pos,t3,by="position")%>%distinct()
      }
      if(filt%in%c("nofilter")){
        t5<-t3
      }
      return(t5)
    }

    gene_list<-sort(unique(cellSTAAR::genes_biomaRt_all%>%filter(.data$gene_biotype=="protein_coding",.data$chromosome_name==!!chr)%>%pull(.data$hgnc_symbol)))
    gene_list<-gene_list[!gene_list==""]
    if(!is.null(genes_manual)){
      gene_list<-genes_manual
    }
    if(grepl("dist",link_type)){
        col_names<-colnames(raw_mappings_dist)
        gene_col<-col_names[grepl("dist",col_names)]
        mappings_cCRE_V3<-bp_level_mappings_dist(raw_mappings_dist%>%filter(!!as.symbol(gene_col)%in%gene_list),filt="nofilter")
      if(nrow(mappings_cCRE_V3)>0){
        index<-logical(length=nrow(mappings_cCRE_V3))
        zzz<-0
        for(pos_check in mappings_cCRE_V3$position){
          zzz<-zzz+1
          index[zzz]<-any(raw_mappings_dist$start<=pos_check & pos_check<=raw_mappings_dist$end)
        }
        if(mean(index)!=1){print("Dist: Some positions do not belong");1+"e"}else{print("Dist: All positions belong")}
      }

    }
    if(grepl("SCREEN_link",link_type)){
      mappings_cCRE_V3<-bp_level_mappings(raw_mappings_SCREEN%>%filter(.data$gene%in%gene_list),filt="nofilter")

      if(nrow(mappings_cCRE_V3)>0){
        index<-logical(length=nrow(mappings_cCRE_V3))
        zzz<-0
        for(pos_check in mappings_cCRE_V3$position){
          zzz<-zzz+1
          index[zzz]<-any(raw_mappings_SCREEN$start<=pos_check & pos_check<=raw_mappings_SCREEN$end)
        }
        if(mean(index)!=1){print("SCREEN: Some positions do not belong");1+"e"}else{print("SCREEN: All positions belong")}
      }
    }
    if(link_type=="EpiMap_link"){
      mappings_cCRE_V3<-bp_level_mappings_EpiMap_link(raw_mappings_EpiMap%>%filter(.data$EpiMap_gene%in%gene_list),filt="nofilter")
        if(nrow(mappings_cCRE_V3)>0){
          index<-logical(length=nrow(mappings_cCRE_V3))
          zzz<-0
          for(pos_check in mappings_cCRE_V3$position){
            zzz<-zzz+1
            index[zzz]<-any(raw_mappings_EpiMap$start<=pos_check & pos_check<=raw_mappings_EpiMap$end)
          }
          if(mean(index)!=1){print("EpiMap: Some positions do not belong");1+"e"}else{print("EpiMap: All positions belong")}
        }
    }
    if(link_type=="ABC_link"){
            mappings_cCRE_V3<-bp_level_mappings_ABC_link(raw_mappings_ABC%>%
                                                           filter(.data$ABC_gene%in%gene_list)
                                                         ,filt="nofilter")
      if(nrow(mappings_cCRE_V3)>0){
        index<-logical(length=nrow(mappings_cCRE_V3))
        zzz<-0
        for(pos_check in mappings_cCRE_V3$position){
          zzz<-zzz+1
          index[zzz]<-any(raw_mappings_ABC$start<=pos_check & pos_check<=raw_mappings_ABC$end)
        }
        if(mean(index)!=1){print("ABC: Some positions do not belong");1+"e"}else{print("ABC: All positions belong")}
      }
    }

    map_obj<-mappings_cCRE_V3
    obj_name<-"mappings_cCRE_V3"
    if(nrow(map_obj)>0){
        if(link_type%in%c("EpiMap_link")){
          cCRE_cols<-which(grepl("EpiMap_gene",colnames(map_obj)))
        }else{
          if(link_type%in%c("ABC_link")){
            cCRE_cols<-which(grepl("ABC_gene",colnames(map_obj)))
          }else{
            cCRE_cols<-which(grepl("cCRE_gene",colnames(map_obj)))
          }

        }

        class_column<-which(colnames(map_obj)=="classification1")
        temp_fun<-function(x,gene){
          return(gene%in%x)
        }
        check_for_matches_bygene<-function(chunk,map_obj,variant_class){
          temp<-matrix(NA,nrow=nrow(map_obj),ncol=length(chunk))
          j<-0
          for(gene in chunk){
            eQTL_subset<-eQTL_matched_file%>%filter(.data$gene_symbol==gene)
            j<-j+1
            temp[,j]<-(map_obj[,class_column]==variant_class &apply(map_obj[,cCRE_cols],MARGIN=1,FUN=temp_fun,gene=gene) & map_obj$cCRE_accession%in%eQTL_subset$cCRE_accession)
          }
          names(temp)<-chunk
          temp<-Matrix(temp,sparse=TRUE)
          return(temp)
        }

        size=50
        chunks<-split(gene_list,ceiling(seq_along(gene_list)/size))
        nchunks<-length(chunks)
        for(z in element_class){
          cl<-parallel::makeForkCluster(ncores)
          registerDoParallel(cl)
          print(paste0("Starting variant type ",z))
          browser()
          print(system.time({
            temp<-pblapply(chunks,FUN=check_for_matches_bygene,map_obj=map_obj,variant_class=z)}))
          parallel::stopCluster(cl)
          temp2<-Matrix(do.call(cbind,temp),sparse=TRUE)
          colnames(temp2)<-gene_list
          rownames(temp2)<-map_obj$position

          out_name<-paste0("variant_mappings_cCRE_V3_",link_type,"_",z,"_","eQTL","_chr",chr)
          assign(eval(out_name),temp2)

          save(list=eval(out_name,envir=environment()),file=paste0(out_wd,"/",out_name,".RData"))
          rm(temp,temp2,mappings_cCRE_V3)
          gc()
        }
        gc()
      }
  } # end link_types loop
seqClose(genofile)
}
