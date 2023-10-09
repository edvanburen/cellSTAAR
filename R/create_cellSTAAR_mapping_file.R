##' @export create_cellSTAAR_mapping_file
create_cellSTAAR_mapping_file<-function(gds.path
                                        ,ct_name
                                        ,num_ct_samples=NULL
                                        ,chr
                                        ,link_types_to_run
                                        ,file_path
                                        ,element_class
                                        ,out_wd
                                        ,genes_manual
                                        ,sc_cutoff=.8){

  genofile <- seqOpen(gds.path)

  filter <- seqGetData(genofile, "annotation/filter")
  AVGDP <- seqGetData(genofile, "annotation/info/AVGDP")
  SNVlist <- filter == "PASS" & AVGDP > 10 & isSNV(genofile)

  variant_pos<-seqGetData(genofile,"position")[SNVlist]%>%enframe()%>%dplyr::rename(position=value)%>%dplyr::select(position)%>%distinct()


  process_bw<-function(path,samp_num=NULL
                       ,ct,chr_filter){
    if(!is.null(samp_num)){
      obj<-as.data.frame(import.bw(paste0(file_path,ct_name,"_",samp_num,".bw")))
    }else{
      obj<-as.data.frame(import.bw(paste0(file_path,ct_name,".bw")))
    }
    obj$seqnames<-as.character(obj$seqnames)
    obj<-obj%>%filter(seqnames==chr_filter)%>%distinct(seqnames,start,end,.keep_all = TRUE)

    obj$position<-0
    t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(seqnames,start,end)%>%mutate(across(position,~.+0:(n() - 1)))%>%ungroup()
    t1$position<-t1$position+t1$start
    t2<-t1%>%dplyr::select(position,score)
    if(!is.null(samp_num)){
      colnames(t2)[2]<-paste0("score_",ct,"_",samp_num)
    }else{
      colnames(t2)[2]<-paste0("score_",ct)
    }

    t2<-right_join(t2,variant_pos,by="position")
    miss<-is.na(t2[,2])
    #t2[miss,2]<-min(t2[!miss,2])/2
    t2[miss,2]<-0
    return(t2)
  }
  process_bed<-function(path,samp_num=NULL,ct,chr_filter){
    if(!is.null(samp_num)){
      obj<-read_delim(paste0(file_path,ct_name,"_",samp_num,".bed.gz"),delim="\t",col_names=FALSE)
    }else{
      obj<-read_delim(paste0(file_path,ct_name,".bed.gz"),delim="\t",col_names=FALSE)
    }

    colnames(obj)<-c("seqnames","start","end","name","qval_int_score","strand","FC","neg_log10pval","neg_log10qval","summit_pos","peak")
    obj<-obj%>%mutate(peak=ct)%>%filter(seqnames==paste0("chr",chr))
    return(obj)
  }

  ### Read in single-cell ATAC-seq data
  if(!is.null(num_ct_samples)){
    ##### If cell type has multiple samples
    ## Peaks
    j<-0
    peak_file<-tibble()
    for(samp_num in 1:num_ct_samples){
      j<-j+1
      peak_file<-bind_rows(peak_file,process_bed(path=file_path
                                                 ,samp_num=samp_num
                                                 ,ct=ct_name
                                                 ,chr_filter = paste0("chr",chr)))
    }
    assign(paste0(ct_name,"_peak"),peak_file%>%
             arrange(seqnames,start,end,FC)%>%distinct(seqnames,start,end,.keep_all = TRUE))

    # Flanking Regions
    j<-0
    region_file<-tibble()
    for(samp_num in 1:num_ct_samples){
      j<-j+1
      if(j==1){
        region_file<-process_bw(path=file_path,samp_num=j
                                ,ct=ct_name,chr_filter = paste0("chr",chr))
      }else{
        new<-process_bw(path=file_path,samp_num=samp_num
                        ,ct=ct_name,chr_filter = paste0("chr",chr))
        region_file<-inner_join(region_file,new,by="position")
      }
    }
    region_file<-region_file%>%mutate(score=rowMeans(across(contains(ct_name))))%>%
      dplyr::select(position,score)%>%ungroup()
    colnames(region_file)[colnames(region_file)=="score"]<-paste0("score_",ct_name)
    assign(ct_name,region_file)
  }else{
    # only 1 sample for the given cell type
    assign(paste0(ct_name,"_peak")
           ,process_bed(path=file_path,ct=ct_name
                ,chr_filter = paste0("chr",chr))%>%
             arrange(seqnames,start,end,FC)%>%distinct(seqnames,start,end,.keep_all = TRUE))

    assign(ct_name,process_bw(path=file_path,ct=ct_name,chr_filter = paste0("chr",chr)))
  }
  #browser()
  for(link_type in link_types_to_run){
    if(link_type=="cCRE_V3_SCREEN_link_eQTL_by_ct"){
      raw_mappings_SCREEN<-agnostic_dnase_summary_V3_eQTL%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,cCRE_accession,gene,.keep_all = TRUE)%>%filter(gene!="")
    }
    browser()
    if(link_type=="cCRE_V3_SCREEN_link_noneQTL_by_ct"){
      raw_mappings_SCREEN<-agnostic_dnase_summary_V3_noneQTL%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,cCRE_accession,gene,.keep_all = TRUE)%>%filter(gene!="")
    }
    if(grepl("dist_0_4000",link_type)){
      raw_mappings_dist<-raw_mappings_cCRE_V3_dist_0_4000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,cCRE_accession,gene_dist_0_4000,.keep_all = TRUE)
    }
    if(grepl("dist_0_1",link_type)){
      raw_mappings_dist<-raw_mappings_cCRE_V3_dist_0_1%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,cCRE_accession,gene_dist_0_1,.keep_all = TRUE)
    }
    if(grepl("dist_1_50000",link_type)){
      raw_mappings_dist<-raw_mappings_cCRE_V3_dist_1_50000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,cCRE_accession,gene_dist_1_50000,.keep_all = TRUE)
    }
    if(grepl("dist_50000_100000",link_type)){
      raw_mappings_dist<-raw_mappings_cCRE_V3_dist_50000_100000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,cCRE_accession,gene_dist_50000_100000,.keep_all = TRUE)
    }
    if(grepl("dist_100000_150000",link_type)){
      raw_mappings_dist<-raw_mappings_cCRE_V3_dist_100000_150000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,cCRE_accession,gene_dist_100000_150000,.keep_all = TRUE)
    }
    if(grepl("dist_150000_200000",link_type)){
      raw_mappings_dist<-raw_mappings_cCRE_V3_dist_150000_200000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,cCRE_accession,gene_dist_150000_200000,.keep_all = TRUE)
    }
    if(grepl("dist_200000_250000",link_type)){
      raw_mappings_dist<-raw_mappings_cCRE_V3_dist_200000_250000%>%filter(chr==paste0("chr",!!chr))
      raw_mappings_dist<-raw_mappings_dist%>%distinct(chr,start,end,cCRE_accession,gene_dist_200000_250000,.keep_all = TRUE)
    }
    if(link_type=="cCRE_V3_EpiMap_link_by_ct"){
      raw_mappings_EpiMap<-raw_mappings_cCRE_V3_EpiMap_link_all_50%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,cCRE_accession,EpiMap_gene,.keep_all = TRUE)
    }
    if(link_type=="cCRE_V3_ABC_link_by_ct"){
      raw_mappings_ABC<-raw_mappings_cCRE_V3_ABC_link_all_50%>%filter(chr==paste0("chr",!!chr))%>%distinct(chr,start,end,cCRE_accession,ABC_gene,.keep_all = TRUE)
    }

    col_names<-colnames(get(paste0(ct_name)))
    col_name<-col_names[grepl("score",col_names)]
    quan<-quantile(get(paste0(ct_name))[,col_name],sc_cutoff,na.rm=TRUE)
    temp_obj<-get(paste0(ct_name))
    #browser()
    ct_CATlas_pos_bw<-get(paste0(ct_name))%>%filter(as.logical(temp_obj[,col_name]>=quan &temp_obj[,col_name]>0))

    all<-get(paste0(ct_name,"_peak"))%>%arrange(chr,start,end)%>%group_by(seqnames,start,end)%>%mutate(num_ct=n(),peak_ct=paste(peak,collapse=","))%>%dplyr::select(-peak)%>%distinct()%>%arrange(seqnames,start,end)%>%ungroup()
    all$width<-all$end-all$start+1
    all$position<-0
    t1<-all%>%dplyr::slice(rep(1:nrow(all),all$width))%>%group_by(seqnames,start,end)%>%mutate(across(position,~.+0:(n() - 1)))%>%ungroup()
    t1$position<-t1$position+t1$start
    ct_CATlas_pos_peak<-t1%>%dplyr::select(position)
    ct_CATlas_pos<-bind_rows(ct_CATlas_pos_peak,ct_CATlas_pos_bw)%>%distinct()

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

      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(cCRE_accession,gene)%>%mutate(across(position,~.+0:(n() - 1)))%>%ungroup()
      t1$position<-t1$position+t1$start

      ##If a position is in multiple regions we really need to remove one
      #temp_to_keep<-t1%>%arrange(position,cCRE_accession)%>%distinct(position,cCRE_accession)%>%group_by(position)%>%dplyr::slice(1)
      #t1<-right_join(t1,temp_to_keep,by=c("cCRE_accession","position"))
      t2<-t1%>%arrange(position)%>%dplyr::select(position,gene,classification1)%>%group_by(position)%>%mutate(genes=paste(gene,collapse=","))%>%dplyr::select(position,genes,classification1)%>%arrange(position)%>%distinct(position,genes,classification1)
      t3<-t2%>%bind_cols(split_into_multiple(column=.$genes,pattern=",",into_prefix="cCRE_gene"))
      colnames(t3)[2]<-paste0("genes_cCRE")

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
    bp_level_mappings_dist<-function(obj,filt,other_cts=FALSE){
      #This is protection in case a group tibble is passed in
      # Running this code on a grouped tibble will cause silent errors!!
      obj<-ungroup(obj)
      if(nrow(obj)==0){return(obj)}
      obj$position<-0
      dist_val<-gsub("cCRE_V3_dist_","",gsub("_by_ct","",link_type))
      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(cCRE_accession,!!as.symbol(paste0("gene_dist_",dist_val)))%>%mutate(across(position,~.+0:(n() - 1)))%>%ungroup()
      t1$position<-t1$position+t1$start

      t2<-t1%>%arrange(position)%>%dplyr::select(position,!!as.symbol(paste0("gene_dist_",dist_val)),classification1)%>%distinct()
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

      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(cCRE_accession,EpiMap_gene)%>%mutate(across(position,~.+0:(n() - 1)))%>%ungroup()

      t1$position<-t1$position+t1$start

      ##If a position is in multiple regions we really need to remove one
      #temp_to_keep<-t1%>%arrange(position,cCRE_accession)%>%distinct(position,cCRE_accession)%>%group_by(position)%>%dplyr::slice(1)
      #t1<-right_join(t1,temp_to_keep,by=c("cCRE_accession","position"))
      t2<-t1%>%arrange(position)%>%dplyr::select(position,EpiMap_gene,classification1)%>%group_by(position)%>%mutate(genes=paste(EpiMap_gene,collapse=","))%>%dplyr::select(position,genes,classification1)%>%arrange(position)%>%distinct(position,genes,classification1)
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

      t1<-obj%>%dplyr::slice(rep(1:nrow(obj),obj$width))%>%group_by(cCRE_accession,ABC_gene)%>%mutate(across(position,~.+0:(n() - 1)))%>%ungroup()
      t1$position<-t1$position+t1$start

      t2<-t1%>%arrange(position)%>%dplyr::select(position,ABC_gene,classification1)%>%group_by(position)%>%mutate(genes=paste(ABC_gene,collapse=","))%>%dplyr::select(position,genes,classification1)%>%arrange(position)%>%distinct(position,genes,classification1)
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

    gene_list<-sort(unique(genes_biomaRt_all%>%filter(gene_biotype=="protein_coding",chromosome_name==chr)%>%pull(hgnc_symbol)))
    gene_list<-gene_list[!gene_list==""]
    if(!is.null(genes_manual)){
      gene_list<-genes_manual
    }
    if(grepl("dist",link_type)){
        col_names<-colnames(raw_mappings_dist)
        gene_col<-col_names[grepl("dist",col_names)]

        mappings_cCRE_V3<-bp_level_mappings_dist(raw_mappings_dist%>%filter(!!as.symbol(gene_col)%in%gene_list),filt="CATlas")
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
    if(grepl("cCRE_V3_SCREEN_link",link_type)){
      mappings_cCRE_V3<-bp_level_mappings(raw_mappings_SCREEN%>%
                                                filter(gene%in%gene_list),filt="CATlas")
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
    if(link_type=="cCRE_V3_EpiMap_link_by_ct"){
          mappings_cCRE_V3<-bp_level_mappings_EpiMap_link(raw_mappings_EpiMap%>%
                                                            filter(EpiMap_gene%in%gene_list)
                                                          ,filt="CATlas")
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
    if(link_type=="cCRE_V3_ABC_link_by_ct"){
          mappings_cCRE_V3<-bp_level_mappings_ABC_link(raw_mappings_ABC%>%
                                                         filter(ABC_gene%in%gene_list)
                                                       ,filt="CATlas")

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
        if(link_type%in%c("cCRE_V3_EpiMap_link_by_ct")){
          cCRE_cols<-which(grepl("EpiMap_gene",colnames(map_obj)))
        }else{
          if(link_type%in%c("cCRE_V3_ABC_link_by_ct")){
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
            j<-j+1
            temp[,j]<-map_obj[,class_column]==variant_class &apply(map_obj[,cCRE_cols],MARGIN=1,FUN=temp_fun,gene=gene)
          }
          names(temp)<-chunk
          temp<-Matrix(temp,sparse=TRUE)
          return(temp)
        }

        size=50
        chunks<-split(gene_list,ceiling(seq_along(gene_list)/size))
        nchunks<-length(chunks)
        for(z in element_class){
          #for(z in c("dELS","pELS")){
          print(paste0("Starting variant type ",z))
          print(system.time({
            temp<-pblapply(chunks,FUN=check_for_matches_bygene,map_obj=map_obj,variant_class=z)}))
          temp2<-Matrix(do.call(cbind,temp),sparse=TRUE)
          colnames(temp2)<-gene_list
          rownames(temp2)<-map_obj$position

          # for some chromosomes
          # this line causes an error when using apply
          # since we are looping over columns
          #
          #temp_summary<-apply(temp2,MARGIN=2,FUN=sum)
          #browser()
          if(grepl("dist",link_type)){
            dist_val<-gsub("cCRE_V3_dist_","",link_type)
            dist_val<-gsub("_by_ct","",dist_val)
            out1_name<-paste0("new2variant_",obj_name,"_new_by_ct_",z,"_",ct_name
                              ,"_dist_",dist_val,"_filter_"
                              ,"CATlas","_",sc_cutoff,"_chr",chr)

          }
          if(grepl("SCREEN",link_type)){
            lt_name<-gsub("_by_ct","",gsub("cCRE_V3_","",link_type))
            out1_name<-paste0("new2variant_",obj_name,"_new_by_ct_",z,"_",ct_name,"_"
                              ,lt_name,"_filter_"
                              ,"CATlas","_",sc_cutoff,"_chr",chr)
          }
          if(link_type=="cCRE_V3_EpiMap_link_by_ct"){
            out1_name<-paste0("new2variant_",obj_name,"_new_by_ct_",z,"_",ct_name,"_EpiMap_link_filter_"
                              ,"CATlas","_",sc_cutoff,"_chr",chr)
          }
          if(link_type=="cCRE_V3_ABC_link_by_ct"){
            out1_name<-paste0("new2variant_",obj_name,"_new_by_ct_",z,"_",ct_name,"_ABC_link_filter_"
                              ,"CATlas","_",sc_cutoff,"_chr",chr)
          }
          assign(eval(out1_name),temp2)
          save(list=eval(out1_name,env=environment()),file=paste0(out_wd,"chr",chr,"/",out1_name,".RData"))
          rm(temp,temp2)
          gc()
        }
        #rm(mappings_cCRE)
        gc()
      }



  }# end link_types loop
seqClose(genofile)
}
