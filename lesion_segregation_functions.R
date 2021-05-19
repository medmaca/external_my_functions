##SET OF FUNCTIONS DESIGNED FOR THE "Lesion_segregation_mutation_summaries.R" SCRIPT AND THE ANALYSIS

#Write a vcf file for reading into MutationalPatterns
write.vcf=function(details,vcf_path,select_vector=NULL,vcf_header_path="~/Documents/vcfHeader.txt") {
  if(class(details)=="character") {
    mat=as.data.frame(stringr::str_split(details,pattern="-",simplify = T),stringsAsFactors=F)
    colnames(mat)<-c("Chrom","Pos","Ref","Alt")
    vcf=create_vcf_files(mat=mat,select_vector=select_vector)
  } else {
    vcf=create_vcf_files(mat=details,select_vector=select_vector)
  }
  write.table(vcf,sep = "\t", quote = FALSE,file=paste0(vcf_path,".temp"),row.names = F)
  system(paste0("cat ",vcf_header_path," ",vcf_path,".temp > ",vcf_path))
  system(paste0("rm ",vcf_path,".temp"))
}

#This function is required in the filtering function
get_ancestor_node=function(node,tree,degree=1){ #to get the 1st degree ancestor (i.e. the direct parent) use degree=1.  Use higher degrees to go back several generations.
  curr<-node
  for(i in 1:degree){
    curr=tree$edge[which(tree$edge[,2]==curr),1]
    if(curr==(1+length(tree$tip.label))) {stop(return(curr))}
  }
  return(curr)
}

#Function for Nick's alternative rho estimation approach (faster)
loglik=function(par,nmuts,depth){
  idx=which(!is.na(nmuts/depth))
  sum(VGAM::dbetabinom(x =nmuts[idx],size = depth[idx],prob = par[2],rho = par[1],log = T))
}

findrho=function(NV_vec,NR_vec){
  pseudo=1e-6
  optres=optim(par=c(0.1,0.1),loglik, gr = NULL,method="L-BFGS-B",lower=c(pseudo,pseudo),upper=c(0.8,0.9),control=list(fnscale=-1),nmuts=NV_vec,depth=NR_vec)
  #c(rho=optres$par[1],p=optres$par[2])
  optres$par[1]
}

#Define the "get_node_types" function required for following the lesion journey in the case of PVVs
#It uses the mutation dataframe ("mut_df") to work out whether daughter branches of a node are (a) a mutant allele (b) wild-type or (c) mixed
get_node_types=function(lesion_children,mut_df,tree) {
  types=sapply(lesion_children, function(node) {
    nodes=c(node,get_all_node_children(node,tree))
    if(all(!mut_df$neg_test[mut_df$clades%in%nodes])){
      return("pure_positive")
    } else if(all(!mut_df$pos_test[mut_df$clades%in%nodes])) {
      return("pure_negative")
    } else {
      return("mixed")
    }
  })
  return(types)
}

#The function to assess all the mutations for whether they are phylogeny breaking
create_PVV_filter_table=function(mutations_to_test,details,tree,matrices,look_back=3,remove_duplicates=F,duplicate_samples=NULL,MC_CORES=1) {
  require(dplyr)
  require(parallel)
  print(paste("Assessing",length(mutations_to_test),"mutations for whether they are truly phylogeny breaking"))
  if(look_back=="all"){
    all_clades=unique(tree$edge[,2])
    all_clade_nodes=lapply(all_clades,function(node) c(node,get_all_node_children(node,tree)))
  }
  filter_output=mclapply(mutations_to_test,function(mut) {
    if(which(mutations_to_test==mut)%%1000==0){print(which(mutations_to_test==mut))}
    allocated_node=details$node[details$mut_ref==mut]
    
    #Get counts of all individual clades within the allocated node i.e. those expected to be positive
    positive_clades=get_all_node_children(allocated_node,tree)
    if(remove_duplicates) {positive_clades=positive_clades[!positive_clades%in%duplicate_samples]} #don't assess duplicate samples as individual samples
    positive_clade_nodes=lapply(positive_clades,function(node) c(node,get_all_node_children(node,tree)))
    positive_clade_samples=lapply(positive_clade_nodes,function(nodes) return(tree$tip.label[nodes[nodes%in%1:length(tree$tip.label)]]))
    NV_pos=unlist(lapply(positive_clade_samples,function(samples) sum(matrices$NV[mut,samples])))
    NR_pos=unlist(lapply(positive_clade_samples,function(samples) sum(matrices$NR[mut,samples])))
    
    #If we previously removed duplicates, need to rederive the "positive_clade_nodes" including duplicates for filtering out ALL positive clades in the next section
    if(remove_duplicates) { 
      positive_clades=get_all_node_children(allocated_node,tree)
      positive_clade_nodes=lapply(positive_clades,function(node) c(node,get_all_node_children(node,tree)))
    }
    
    #Get counts of nearby individual clades that don't include the allocated node, i.e. those expected to be negative
    if(look_back=="all"){
      select=lapply(all_clade_nodes,function(nodes)if(!any(nodes%in%positive_clade_nodes)){TRUE}else{FALSE})
      negative_clades=all_clades[unlist(select)]
      if(remove_duplicates) {negative_clades=negative_clades[!negative_clades%in%duplicate_samples]}
      negative_clade_nodes=lapply(negative_clades,function(node) c(node,get_all_node_children(node,tree)))
      negative_clade_samples=lapply(negative_clade_nodes,function(nodes) return(tree$tip.label[nodes[nodes%in%1:length(tree$tip.label)]]))
      NV_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NV[mut,samples])))
      NR_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NR[mut,samples])))
    } else {
      #Assess only "nearby" to the allocated node (up to the degree specified by look_back argument)
      ancestral_node=get_ancestor_node(allocated_node,tree,degree=look_back)
      all_clades=c(ancestral_node,get_all_node_children(ancestral_node,tree))
      all_clade_nodes=lapply(all_clades,function(node) c(node,get_all_node_children(node,tree)))
      
      select=lapply(all_clade_nodes,function(nodes)if(!any(nodes%in%positive_clade_nodes)){TRUE}else{FALSE})
      negative_clades=all_clades[unlist(select)]
      if(remove_duplicates) {negative_clades=negative_clades[!negative_clades%in%duplicate_samples]}
      negative_clade_nodes=lapply(negative_clades,function(node) c(node,get_all_node_children(node,tree)))
      negative_clade_samples=lapply(negative_clade_nodes,function(nodes) return(tree$tip.label[nodes[nodes%in%1:length(tree$tip.label)]]))
      NV_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NV[mut,samples])))
      NR_neg=unlist(lapply(negative_clade_samples,function(samples) sum(matrices$NR[mut,samples])))
    }
    
    #Apply beta-binomial filter to these counts to get sense of overdispersion
    NR_pos[NR_pos==0]<-1; NR_neg[NR_neg==0]<-1 #Set sites with 0 depth to a depth of 1
    pos_rho=estimateRho_gridml(NV_vec=NV_pos,NR_vec=NR_pos)
    neg_rho=estimateRho_gridml(NV_vec=NV_neg,NR_vec=NR_neg)
    
    #Additional filters to check not just for overdispersion, but that at least one clade unexpectedly truly negative or positive
    pos_test=any(NV_pos==0 & NR_pos>=10) #is there at least one "WT" sub-clade with zero variant reads with a depth of ≥ 10
    neg_test=any((NV_neg/NR_neg)>=0.3 & NR_neg >=8) #is there at least one anticipated "negative" clade, that in fact has a VAF>0.3 with a depth ≥8
    
    mut_params=data.frame(mut=mut,
                          node=allocated_node,
                          pos_rho=pos_rho,
                          neg_rho=neg_rho,
                          pos_test=pos_test,
                          neg_test=neg_test,
                          pval=details$pval[which(details$mut_ref==mut)])
    return(mut_params) 
  },mc.cores = MC_CORES)
  
  filter_output_df=dplyr::bind_rows(filter_output)
  return(filter_output_df)
}

#As above, but can incorporate two alternative mutant alleles, therefore suitable for MAVs
get_MAV_node_types=function(lesion_children,mut_df,tree) {
  types=sapply(lesion_children, function(node) {
    nodes=c(node,get_all_node_children(node,tree))
    if(all(!mut_df$neg_test[mut_df$clades%in%nodes])&any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])&!any(mut_df$mut2_pos_test[mut_df$clades%in%nodes])){
      return("pure_mut1")
    } else if((sum(!mut_df$neg_test[mut_df$clades%in%nodes])/length(!mut_df$neg_test[mut_df$clades%in%nodes])>0.98)&any(mut_df$mut2_pos_test[mut_df$clades%in%nodes])&!any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])){
      return("pure_mut2")
    }else if(all(!mut_df$mut1_pos_test[mut_df$clades%in%nodes]) & all(!mut_df$mut2_pos_test[mut_df$clades%in%nodes])) {
      return("pure_negative")
    } else {
      return("mixed")
    }
  })
  return(types)
}

#Function to reclassify nearby mutations that are on the same branch as multi-nucleotide variants
#Note - this function ASSUMES that the mutations are phased, and this should be checked
#Needs a reference genome file
reclassify_MNVs=function(COMB_mats,region_size=2,genomeFile) {
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  library("stringr")
  
  details=COMB_mats$mat
  NV=COMB_mats$NV
  NR=COMB_mats$NR
  
  Chroms=c(1:22,"X","Y")
  out_list_by_chrom=lapply(Chroms,function(Chrom) {
    details_by_chrom=details[details$Chrom==as.character(Chrom) & details$Mut_type=="SNV",]
    out_list=lapply(1:nrow(details_by_chrom),function(i) {
      Pos=as.numeric(details_by_chrom$Pos[i])
      node=details_by_chrom$node[i]
      near_muts=which(as.numeric(details_by_chrom$Pos)>Pos&
                        as.numeric(details_by_chrom$Pos)<=(Pos+region_size)&
                        details_by_chrom$node==node)
      near_muts<-near_muts[near_muts!=i]
      if(length(near_muts)==0) {
        return(NA)
      } else {
        return(c(details_by_chrom$mut_ref[i],details_by_chrom$mut_ref[near_muts]))
      }
    })
    out_list[sapply(out_list,function(x) is.na(x[1]))]<-NULL
    return(out_list)
  })
  
  out_list_by_chrom=lapply(out_list_by_chrom, function(list) {
    i=1
    while(i<length(list)){
      muts<-list[[i]]
      if(any(list[[i+1]]%in%muts)) {list[[i]]<-unique(c(list[[i]],list[[i+1]]));list[[i+1]]<-NULL} else {i<-i+1}
    }
    return(list)
  })
  
  #Combine the chromosomes into one list
  out_list=unlist(out_list_by_chrom,recursive=F)
  
  print(paste("There are",length(out_list),"SNV pairs that will be reclassified as MNVs"))
  
  # #Trinucleotide substitutions will occur more than once - therefore need to replace these (not yet done)
  if(length(out_list)>0) {
    new_df=lapply(out_list,function(x) {
      chrom=str_split(x[1],pattern="-",simplify=T)[,1]
      pos=as.numeric(str_split(x,pattern="-",simplify=T)[,2])
      ref=str_split(x,pattern="-",simplify=T)[,3]
      alt=str_split(x,pattern="-",simplify=T)[,4]
      node=details$node[details$mut_ref==x[1]]
      pval=details$pval[details$mut_ref%in%x]
      if(length(pos)==2 & pos[2]==(pos[1]+1)){
        new_ref=paste(ref,collapse="")
        new_alt=paste(alt,collapse="")
        return(data.frame(mut_ref=paste(chrom,pos[1],new_ref,new_alt,sep="-"),Chrom=chrom,Pos=pos[1],Ref=new_ref,Alt=new_alt,Mut_type="MNV",node=node,pval=mean(pval)))
      } else {
        new_ref=as.vector(scanFa(genomeFile, GRanges(chrom, IRanges(min(pos),max(pos)))))
        new_ref_vec=as.vector(str_split(new_ref,"",simplify=T))
        names(new_ref_vec)=seq(min(pos),max(pos),by=1)
        new_alt_vec=new_ref_vec
        for(i in 1:length(pos)){
          new_alt_vec[as.character(pos[i])]<-alt[i]
        }
        new_alt=paste(new_alt_vec,collapse="")
        return(data.frame(mut_ref=paste(chrom,pos[1],new_ref,new_alt,sep="-"),Chrom=chrom,Pos=pos[1],Ref=new_ref,Alt=new_alt,Mut_type="MNV",node=node,pval=mean(pval)))
      }
    })
    new_df=dplyr::bind_rows(new_df)
    
    replaced_SNVs=unlist(out_list)
    SNVs_for_counts=sapply(out_list,function(x) x[1])
    
    new_NV=NV[SNVs_for_counts,]
    new_NR=NR[SNVs_for_counts,]
    
    rownames(new_NV)=rownames(new_NR)<-new_df$mut_ref
    
    details_new<-details[!details$mut_ref%in%replaced_SNVs,]
    details_new<-rbind(details_new[,colnames(new_df)],new_df)
    
    return(list(mat=details_new,NV=rbind(NV[!details$mut_ref%in%replaced_SNVs,],new_NV),NR=rbind(NR[!details$mut_ref%in%replaced_SNVs,],new_NR)))
  } else {
    return(COMB_mats)
  }
}

#Function to find potential multi-allelic variants that may be caused by persistent DNA lesions
#Assesses for overlapping positions of reference
get_multi_allelic_variant_list=function(details,SNV_only=F) {
  Chroms=c(1:22,"X","Y")
  details$Ref=as.character(details$Ref)
  details$Alt=as.character(details$Alt)
  if(SNV_only){
    if(!"Chrom_pos"%in%colnames(details)){
      details$Chrom_pos<-paste(details$Chrom,details$Pos,sep = "-")
    }
    duplicates=details$Chrom_pos[duplicated(details$Chrom_pos)]
    dups_list<-lapply(duplicates,function(Chrom_pos) {return(details$mut_ref[details$Chrom_pos==Chrom_pos])})
    return(dups_list)
  } else {
    out_list_by_chrom=lapply(Chroms,function(Chrom) {
      if(sum(details$Chrom==Chrom)>1) {
        print(paste("Analysing chromosome",Chrom))
        #Split comparison of mutations by chromosome
        details_by_chrom=details[details$Chrom==as.character(Chrom),]
        details_by_chrom=details_by_chrom[order(details_by_chrom$Pos),]
        dups_list=lapply(1:(nrow(details_by_chrom)-1),function(i) {
          if(i%%1000==0) {print(i)}
          if(nchar(details_by_chrom$Ref[i])==1) {
            Pos=as.numeric(details_by_chrom$Pos[i])
          } else {
            Pos=as.numeric(details_by_chrom$Pos[i]):(as.numeric(details_by_chrom$Pos[i])+nchar(details_by_chrom$Ref[i])-1)
          }
          out=sapply((i+1):min(i+5,nrow(details_by_chrom)),function(j) {
            if(nchar(details_by_chrom$Ref[j])==1) {
              pos=as.numeric(details_by_chrom$Pos[j])
            } else {
              pos=as.numeric(details_by_chrom$Pos[j]):(as.numeric(details_by_chrom$Pos[j])+nchar(details_by_chrom$Ref[j])-1)
            }
            if(length(intersect(pos,Pos))>0) {
              return(T)
            } else {
              return(F)
            }
          })
          if(any(out)) {return(c(details_by_chrom$mut_ref[i],details_by_chrom$mut_ref[(i+1):min(i+5,nrow(details_by_chrom))][out]))} else {return(NA)}
        })
        dups_list[sapply(dups_list,function(x) is.na(x[1]))]<-NULL
        return(dups_list)
      } else {
        return(NULL)
      }
      
    })
    out_list=unlist(out_list_by_chrom,recursive=F)
    return(out_list)
  }
}

#Find the latest possible timing of the acquisition of the lesion
find_PVV_lesion_node=function(mut,allocated_node,pos_test,neg_test,tree,matrices) {
  
  #Define the get_ancestor_node function
  get_ancestor_node=function(node,tree,degree=1){ #to get the 1st degree ancestor (i.e. the direct parent) use degree=1.  Use higher degrees to go back several generations.
    curr<-node
    for(i in 1:degree){
      curr=tree$edge[which(tree$edge[,2]==curr),1]
      if(curr==(1+length(tree$tip.label))) {stop(return(curr))}
    }
    return(curr)
  }
  
  mut_df=create_mut_df(mut=mut,tree=tree,matrices=matrices)
  
  if(pos_test) {
    #if the negative sub-clade is within the allocated node, then "pos_test" will be true and "allocated_node" is the "initial_lesion_node"
    initial_lesion_node<-allocated_node
  } else if (neg_test){
    #If there is a positive clade outside the allocated node, then "neg_test" will be true.  In this case need to find the positive clade.
    #Do this be iteratively going from the allocated node to its ancestral node and looking for the ancestral node that contains ALL positive clades in the tree
    all_clades=unique(tree$edge[,2])
    all_clade_samples=lapply(all_clades,function(node) getTips(node=node,tree=tree))
    
    #Iteratively look back through ancestral nodes to find the one encasing other positive samples
    j=1
    repeat {
      ancestor=get_ancestor_node(allocated_node,tree,degree=j)
      ancestor_tips=getTips(tree,ancestor)
      #Look at clades that don't have any of the samples in "ancestor_tips". If this ancestor is the "initial_lesion_node", none will meet the "pos_test".
      if(!any(mut_df$pos_test[unlist(lapply(all_clade_samples, function(samples) !any(samples %in% ancestor_tips)))])) {
        break #Once this criteria is met, do not need to look back any further
      }
      else if(ancestor==tree$edge[1,1]){
        break
      }
      j=j+1
    }
    initial_lesion_node<-ancestor #The initial lesion node is therefore the most recent "ancestor" from the previous loop
  }
  return(initial_lesion_node)
}

find_MAV_lesion_node=function(node1,node2,tree,Chrom="auto") {
  if(node1==node2 & !Chrom%in%c("X","Y")) {
    stop(return(return(list(initial_lesion_node=NA,Filter="FAIL",Class="FAIL"))))
  } else {
    #get the ancestral nodes
    ancestor_1=tree$edge[tree$edge[,2]==node1,1]
    ancestor_2=tree$edge[tree$edge[,2]==node2,1]
    ancestors=c(ancestor_1,ancestor_2)
    ancestor_heights=sapply(ancestors,function(node) nodeheight(tree,node))
    
    if(ancestor_1==ancestor_2){
      Filter="PASS"
      Class="simple"
      initial_lesion_node=get_ancestor_node(node = node1,tree = tree)
    } else if(ancestor_1%in%get_all_node_children(ancestor_2,tree)|ancestor_2%in%get_all_node_children(ancestor_1,tree)){
      Filter="PASS"
      Class="removed"
      if(which.min(ancestor_heights)==1 & node2%in%get_all_node_children(node1,tree)) {
        initial_lesion_node<-node1
      } else if(which.min(ancestor_heights)==2 & node1%in%get_all_node_children(node2,tree)){
        initial_lesion_node<-node2
      } else {
        initial_lesion_node<-ancestors[which.min(ancestor_heights)]
      }
    } else {
      return(return(list(initial_lesion_node=NA,Filter="FAIL",Class="FAIL")))
    }
  }
  return(list(initial_lesion_node=initial_lesion_node,Filter=Filter,Class=Class))
}

#This function extracts a straight-forward phasing info df from the julia algorithm output
extract_phasing_info=function(list,Ref,Alt) {
  phasing_df=dplyr::bind_rows(list)
  #basects_df=Reduce(rbind,lapply(list,function(list) return(list[[2]])))
  if(is.logical(phasing_df)) {
    stop(return(NA))
  } else if(nrow(phasing_df)==0){
    stop(return(NA))
  }
  phasing_df$mut_base=sapply(strsplit(phasing_df$Mutation_allele,split = "="),function(x) x[2])
  phasing_df$snp_base=sapply(strsplit(phasing_df$SNP_allele,split = "="),function(x) x[2])
  
  SNP_sites=sort(unique(str_split(phasing_df$SNP_allele,pattern = "=",simplify = T)[,1]))
  
  phasing_by_SNP_list=lapply(SNP_sites,function(SNP_site) {
    phasing_df_snp<-phasing_df[grepl(SNP_site,phasing_df$SNP_allele),]
    
    #Summarise mut allele phasing
    if(any(phasing_df_snp$mut_base%in%Alt)) {
      alt_phasing=table(phasing_df_snp$snp_base[phasing_df_snp$mut_base%in%Alt])
      alt_phases_with_base=names(alt_phasing)[which.max(alt_phasing)]
      n_alt_reads_supporting=alt_phasing[alt_phases_with_base]
      n_alt_reads_against=(sum(alt_phasing)-n_alt_reads_supporting)
    } else {
      alt_phasing=NA
      alt_phases_with_base=NA
      n_alt_reads_supporting=0
      n_alt_reads_against=0
    }
    
    #Summarise wt allele phasing
    if(any(phasing_df_snp$mut_base==Ref)){
      ref_phasing=table(phasing_df_snp$snp_base[phasing_df_snp$mut_base==Ref])
      ref_phases_with_base=names(ref_phasing)[which.max(ref_phasing)]
      n_ref_reads_supporting=ref_phasing[ref_phases_with_base]
      n_ref_reads_against=(sum(ref_phasing)-n_ref_reads_supporting)
    } else {
      ref_phasing=NA
      ref_phases_with_base=NA
      n_ref_reads_supporting=0
      n_ref_reads_against=0
    }
    
    return(data.frame(SNP_site=SNP_site,
                      alt_phases_with_base=alt_phases_with_base,
                      n_alt_reads_supporting=n_alt_reads_supporting,
                      n_alt_reads_against=n_alt_reads_against,
                      ref_phases_with_base=ref_phases_with_base,
                      n_ref_reads_supporting=n_ref_reads_supporting,
                      n_ref_reads_against=n_ref_reads_against
    )
    )
  })
  phasing_by_SNP_df=dplyr::bind_rows(phasing_by_SNP_list)
  return(phasing_by_SNP_df)
}

#Function will look in the supplied output_dir to see if phasing output for given sample/Chrom/Pos already exists, if not will run the .jl script. Imports the data.
#Run example: get_phasing_list(samples=positive_samples1,Chrom=Chrom,Pos=Pos,project=project,output_dir = phasing_output_dir,ref_sample_set = Ref_sample_set)

get_phasing_list=function(samples,Chrom,Pos,project,tree=NULL,output_dir,ref_sample_set,distance=1000,force_rerun=F,verbose=F,use_tree=T) {
  wd<-getwd()
  setwd("/lustre/scratch119/realdata/mdt1/team154/ms56/my_programs/Mike_phasing") #Need to be in this directory for the function
  if(is.numeric(project)) {
    phasing_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample));print(paste("Reference sample set chosen as",ref_sample_set))}
      if(!file.exists(phasing_output_file)|file.info(phasing_output_file)$size==0|force_rerun) {
        #This section is to account for the long bam headers in sample PD44579b which interfere with the script
        if(grepl("PD44579b",sample)) {
          #Import all the necessary bams with edited headers
          sapply(c(sample,unlist(strsplit(ref_sample_set,","))),function(bam_sample) {
            new_bam_path=paste0("new_bams/",bam_sample,".sample.dupmarked.bam")
            if(!file.exists(new_bam_path)) {
              print(paste("Importing bam file for",bam_sample,"and replacing header"))
              bam_path=paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",bam_sample,"/",bam_sample,".sample.dupmarked.bam")
              command=paste("julia header_edit.jl",bam_path,"offending_string.txt")
              system(command)
            }
          })
          #Now run using the modified julia script to use these local files
          command=paste("julia DRIVER_phasing_specify_BAM_directory.jl",Chrom,Pos,sample,"/lustre/scratch119/casm/team154pc/ms56/my_programs/Mike_phasing/new_bams",as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        } else {
          command=paste("julia DRIVER_phasing.jl",Chrom,Pos,sample,project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        }
      } else if(verbose) {
        print("Existing phasing files found in specified output directory")
      }
      if(file.exists(phasing_output_file)&file.info(phasing_output_file)$size!=0){
        phasing=read.table(phasing_output_file,header = T,stringsAsFactors = F,colClasses="character")
      } else {
        print("Unable to run phasing script")
        phasing=NA
      }
      return(phasing)
    })
  } else if(is.data.frame(project)) {
    phasing_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample))}
      sample_project=project$project[project$sample==sample]
      if(use_tree){
        set.seed(1)
        ref_sample_set=paste0(sample(x=tree$tip.label[tree$tip.label%in%project$sample[project$project==sample_project]],size=5),collapse=",")
      } else {
        sample_stem=stringr::str_split(sample,pattern = "_",simplify=T)[,1]
        set.seed(1)
        ref_sample_set=paste0(sample(x=project$sample[project$project==sample_project & grepl(sample_stem,project$sample)],size=5),collapse=",")
      }
      if(verbose) {print(paste("Ref sample set chosen as",ref_sample_set))}
      
      if(!file.exists(phasing_output_file)|file.info(phasing_output_file)$size==0|force_rerun) {
        command=paste("julia DRIVER_phasing.jl",Chrom,Pos,sample,sample_project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
        system(command)
      } else if(verbose) {
        print("Existing phasing files found in specified output directory")
      }
      if(file.exists(phasing_output_file)&file.info(phasing_output_file)$size!=0){
        phasing=read.table(phasing_output_file,header = T,stringsAsFactors = F,colClasses="character")
      } else {
        phasing="Unable to run phasing script"
      }
      return(phasing)
    })
  }
  setwd(wd)
  return(phasing_list)
}


get_base_counts_list=function(samples,Chrom,Pos,project,tree=NULL,output_dir,ref_sample_set,distance=1000,force_rerun=F,verbose=F,use_tree=T) {
  wd<-getwd()
  setwd("/lustre/scratch119/realdata/mdt1/team154/ms56/my_programs/Mike_phasing") #Need to be in this directory for the function
  if(is.numeric(project)) {
    basects_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample));print(paste("Reference sample set chosen as",ref_sample_set))}
      if(!file.exists(basects_output_file)|file.info(basects_output_file)$size==0|force_rerun) {
        if(grepl("PD44579b",sample)) {
          #Import all the necessary bams with edited headers
          sapply(c(sample,unlist(strsplit(ref_sample_set,","))),function(bam_sample) {
            new_bam_path=paste0("new_bams/",bam_sample,".sample.dupmarked.bam")
            if(!file.exists(new_bam_path)) {
              print(paste("Importing bam file for",bam_sample,"and replacing header"))
              bam_path=paste0("/nfs/cancer_ref01/nst_links/live/",project,"/",bam_sample,"/",bam_sample,".sample.dupmarked.bam")
              command=paste("julia header_edit.jl",bam_path,"offending_string.txt")
              system(command)
            }
          })
          #Now run using the modified julia script to use these local files
          command=paste("julia DRIVER_phasing_specify_BAM_directory.jl",Chrom,Pos,sample,"/lustre/scratch119/casm/team154pc/ms56/my_programs/Mike_phasing/new_bams",as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        } else {
          command=paste("julia DRIVER_phasing.jl",Chrom,Pos,sample,project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
          system(command) 
        }
      } else if(verbose) {
        print("Existing base counts files found in specified output directory")
      }
      if(file.exists(basects_output_file)&file.info(basects_output_file)$size!=0){
        basects=read.table(basects_output_file,header = T,stringsAsFactors = F)
      } else {
        basects="Unable to run phasing script"
      }
      return(basects)
    })
  } else if(is.data.frame(project)) {
    basects_list=lapply(samples,function(sample) {
      phasing_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_phasing.txt")
      basects_output_file=paste0(output_dir,"/",sample,"_",Chrom,"_",Pos,"_basects.txt")
      if(verbose) {print(paste("Looking in sample",sample))}
      sample_project=project$project[project$sample==sample]
      if(use_tree){
        set.seed(1)
        ref_sample_set=paste0(sample(x=tree$tip.label[tree$tip.label%in%project$sample[project$project==sample_project]],size=5),collapse=",")
      } else {
        sample_stem=stringr::str_split(sample,pattern = "_",simplify=T)[,1]
        set.seed(1)
        ref_sample_set=paste0(sample(x=project$sample[project$project==sample_project & grepl(sample_stem,project$sample)],size=5),collapse=",")
      }
      if(verbose) {print(paste("Ref sample set chosen as",ref_sample_set))}
      if(!file.exists(basects_output_file)|file.info(basects_output_file)$size==0|force_rerun) {
        command=paste("julia DRIVER_phasing.jl",Chrom,Pos,sample,sample_project,as.character(distance),phasing_output_file,basects_output_file,ref_sample_set)
        system(command)
      } else if(verbose) {
        print("Existing base counts files found in specified output directory")
      }
      if(file.exists(basects_output_file)&file.info(basects_output_file)$size!=0){
        basects=read.table(basects_output_file,header = T,stringsAsFactors = F)
      } else {
        basects="Unable to run phasing script"
      }
      return(basects)
    })
  }
  setwd(wd)
  return(basects_list)
}

get_clade_base_counts=function(nodes,tree,Chrom,Pos,project,ref_sample_set,phasing_output_dir,distance=1000,force_rerun=F) {
  clade_base_counts_list=lapply(nodes,function(node) {samples=getTips(tree=tree,node=node);base_counts_list=get_base_counts_list(samples=samples,Chrom=Chrom,Pos=Pos,project=project,tree=tree,output_dir = phasing_output_dir,ref_sample_set = ref_sample_set,distance=distance,force_rerun = force_rerun);return(base_counts_list)})
  aggregated_clade_base_counts_list=lapply(clade_base_counts_list,function(list) {
    chrom_pos_df=list[[1]][,c(2,3)] #Get the co-ordinates of apparent het SNPs from the 1st in the list
    list_mod=lapply(list,function(df) {res<-left_join(chrom_pos_df,df[,-1],by=c("Chr","Pos"));res[is.na(res)]<-0;return(res[,3:7])}) #Now get just the base counts at these sites
    counts_df=Reduce(function(df1,df2) {return(df1+df2)},list_mod) #Aggregate these across a pure subclade
    return(cbind(chrom_pos_df,counts_df))
  })
  return(aggregated_clade_base_counts_list)
}

return_heterozygous_SNPs=function(base_counts_list) {
  pos_het_res=lapply(base_counts_list,function(df) {
    if(nrow(df)==0) {stop(return(NA))}
    het_SNPs=apply(df[,3:7],1,function(x) {
      counts=x; sum_counts=sum(x)
      if(sum_counts==0) {stop(return(NA))}
      het_test=sapply(counts,function(y) return(binom.test(x=y,n=sum_counts)$p.value))
      absent_test=sapply(counts,function(y) return(binom.test(x=y,n=sum_counts,p=0.01)$p.value))
      hom_test=sapply(counts,function(y) return(binom.test(x=y,n=sum_counts,p=0.99)$p.value))
      lik_hom=prod(pmax(hom_test,absent_test))
      lik_het=prod(pmax(het_test,absent_test))
      if(lik_het>lik_hom) {return(T)} else {return(F)}
    })
    return(het_SNPs)
  })
  positions=base_counts_list[[1]]$Pos
  het_positions=positions[Reduce(function(x,y) {x&y},pos_het_res)]
  return(het_positions)
}

check_matching_phasing=function(phasing_info1,phasing_info2,het_positions=NULL) {
  if(any(sapply(list(phasing_info1,phasing_info2),is.logical))) {
    result<-"Unable to confirm phasing"
  } else {
    comb_df=inner_join(phasing_info1,phasing_info2,by="SNP_site")
    comb_df$depth=comb_df$n_alt_reads_supporting.x+comb_df$n_alt_reads_supporting.y+comb_df$n_ref_reads_supporting.x+comb_df$n_ref_reads_supporting.y
    
    #Test if the called SNPs appear real, or if previously tested on the basects data, filter the included SNPs based on this
    if(is.null(het_positions)) {
      true_het1<-comb_df$alt_phases_with_base.x!=comb_df$ref_phases_with_base.x
      true_het2<-comb_df$alt_phases_with_base.y!=comb_df$ref_phases_with_base.y
      het_test=mapply(FUN=function(x,y) {vec=c(x,y);vec<-vec[!is.na(vec)];if(length(vec)==0) {return(T)} else if(all(vec)) {return(T)} else {return(F)}},x=true_het1,y=true_het2)
      
      #If no true het SNPs, stop function; else filter the comb_df
      if(!any(het_test)) {
        stop(return("No SNPs appear to be heterozygous"))
      } else {
        comb_df<-comb_df[het_test,]
      }
      
      #Although heterozygosity is likely after the above test, it is not confirmed. Test for this:
      het_confirmed=apply(comb_df[,c("ref_phases_with_base.x","ref_phases_with_base.y","alt_phases_with_base.x","alt_phases_with_base.y")],1,function(x) length(unique(x[!is.na(x)]))>1)
      if(any(het_confirmed)) {
        comb_df<-comb_df[het_confirmed,]
        het_not_confirmed<-F
      } else {
        het_not_confirmed<-T
      }
      
    } else {
      snp_positions=as.numeric(str_split(pattern=":",comb_df$SNP_site,simplify=T)[,2])
      comb_df<-comb_df[snp_positions%in%het_positions,]
      if(nrow(comb_df)==0) {
        stop(return("No SNPs appear to be heterozygous"))
      }
      het_not_confirmed<-F
    }
    
    
    #Test for either alt matching alt, or ref matching ref for any individual SNP
    matching_res=list(Matching_alt_phasing=comb_df$alt_phases_with_base.x==comb_df$alt_phases_with_base.y,
                      Matching_ref_phasing=comb_df$ref_phases_with_base.x==comb_df$ref_phases_with_base.y,
                      Non_matching_alt_ref_phasing=comb_df$alt_phases_with_base.x!=comb_df$ref_phases_with_base.y,
                      Non_matching_ref_alt_phasing=comb_df$ref_phases_with_base.x!=comb_df$alt_phases_with_base.y)
    
    matching_res=lapply(matching_res, function(vec) {
      names(vec)<-1:length(vec)
      vec_no_NAs<-vec[!is.na(vec)]
        if(length(unique(vec_no_NAs))>1) { #If there is disagreement between different SNPs, retain the highest depth ones only
          vec_no_NAs<-vec_no_NAs[as.character(which(comb_df$depth>median(comb_df$depth)))]
          }
        return(vec_no_NAs)
    })
    
    if(all(sapply(matching_res,function(x) length(x)==0))){
      result<-"Unable to confirm phasing"
    } else if(all(unlist(matching_res))) {
      if(het_not_confirmed) {
        result<-"Same phasing suggested, though SNP heterozygosity not confirmed"
      } else {
        result<-"Same phasing confirmed"
      }
    } else {
      if(het_not_confirmed) {
        result<-"Non-matching suggested, though SNP heterozygosity not confirmed"
      } else {
        result<-"Non-matching phasing confirmed"
      }
    }
  }
  return(result)
}


check_matching_phasing_non_clonal=function(phasing_info1,phasing_info2,het_positions=NULL) {
  if(any(sapply(list(phasing_info1,phasing_info2),is.logical))) {
    result<-"Unable to confirm phasing"
  } else {
    if(is.null(het_positions)){stop(return("Must supply list of confirmed heterozygous SNP positions"))}
    comb_df=inner_join(phasing_info1,phasing_info2,by="SNP_site")
    comb_df$depth=comb_df$n_alt_reads_supporting.x+comb_df$n_alt_reads_supporting.y+comb_df$n_ref_reads_supporting.x+comb_df$n_ref_reads_supporting.y
    snp_positions=as.numeric(str_split(pattern=":",comb_df$SNP_site,simplify=T)[,2])
    comb_df<-comb_df[snp_positions%in%het_positions,]
    if(nrow(comb_df)==0) {
      stop(return("No SNPs appear to be heterozygous"))
    }
    
    #Test for either alt matching alt, or ref matching ref for any individual SNP
    matching_res=list(Matching_alt_phasing=comb_df$alt_phases_with_base.x==comb_df$alt_phases_with_base.y)
    
    matching_res=lapply(matching_res, function(vec) {
      names(vec)<-1:length(vec)
      vec_no_NAs<-vec[!is.na(vec)]
      if(length(unique(vec_no_NAs))>1) { #If there is disagreement between different SNPs, retain the highest depth ones only
        vec_no_NAs<-vec_no_NAs[as.character(which(comb_df$depth>median(comb_df$depth)))]
      }
      return(vec_no_NAs)
    })
    
    if(all(sapply(matching_res,function(x) length(x)==0))){
      result<-"Unable to confirm phasing"
    } else if(all(unlist(matching_res))) {
      result<-"Same phasing confirmed"
    } else {
      result<-"Non-matching phasing confirmed"
    }
  }
  return(result)
}

assess_phasing_non_clonal=function(phasing_summaries,sample_sets,Chrom,Pos,project,tree=NULL,output_dir,ref_sample_set,distance=1000,use_tree=T){
  phase_sum1=phasing_summaries[[1]]
  phase_sum2=phasing_summaries[[2]]
  if(is.logical(phase_sum1)|is.logical(phase_sum2)) {
    stop(return("Unable to phase"))
  } else {
    clade_base_counts_list=lapply(sample_sets,function(samples) {base_counts_list=get_base_counts_list(samples=samples,Chrom=Chrom,Pos=Pos,project=project,tree=tree,output_dir = output_dir,ref_sample_set = ref_sample_set,distance=distance,use_tree=use_tree);return(base_counts_list)})
    aggregated_clade_base_counts_list=lapply(clade_base_counts_list,function(list) {
      chrom_pos_df=list[[1]][,c(2,3)] #Get the co-ordinates of apparent het SNPs from the 1st in the list
      list_mod=lapply(list,function(df) {res<-left_join(chrom_pos_df,df[,-1],by=c("Chr","Pos"));res[is.na(res)]<-0;return(res[,3:7])}) #Now get just the base counts at these sites
      counts_df=Reduce(function(df1,df2) {return(df1+df2)},list_mod) #Aggregate these across a pure subclade
      return(cbind(chrom_pos_df,counts_df))
    })
    het_positions<-return_heterozygous_SNPs(base_counts_list=aggregated_clade_base_counts_list)
    het_positions<-het_positions[het_positions!=Pos] #Exclude the position of the actual mutation
    outcome=check_matching_phasing_non_clonal(phase_sum1,phase_sum2,het_positions = het_positions)
    return(outcome)
  }
}

get_confirmed_heterozygous_SNPs=function(phasing_info1,phasing_info2){
  comb_df=full_join(phasing_info1,phasing_info2,by="SNP_site")
  comb_df$depth=comb_df$n_alt_reads_supporting.x+comb_df$n_alt_reads_supporting.y+comb_df$n_ref_reads_supporting.x+comb_df$n_ref_reads_supporting.y
  
  #Test if the called SNPs appear real
  true_het1<-comb_df$alt_phases_with_base.x!=comb_df$ref_phases_with_base.x
  true_het2<-comb_df$alt_phases_with_base.y!=comb_df$ref_phases_with_base.y
  het_test=mapply(FUN=function(x,y) {vec=c(x,y);vec<-vec[!is.na(vec)];if(length(vec)==0) {return(T)} else if(all(vec)) {return(T)} else {return(F)}},x=true_het1,y=true_het2)
  
  #If no true het SNPs, stop function; else filter the comb_df
  if(!any(het_test)) {
    stop(return("No SNPs appear to be heterozygous"))
  } else {
    comb_df<-comb_df[het_test,]
  }
  
  #Although heterozygosity is likely after the above test, it is not confirmed. Test for this:
  het_confirmed=apply(comb_df[,c("ref_phases_with_base.x","ref_phases_with_base.y","alt_phases_with_base.x","alt_phases_with_base.y")],1,function(x) length(unique(x[!is.na(x)]))>1)
  if(any(het_confirmed)) {
    comb_df<-comb_df[het_confirmed,]
    return(comb_df$SNP_site)
  } else {
    return(NA)
  }
}


check_for_both_alleles_confirming_ref=function(phasing_info) {
  if(is.logical(phasing_info)) {
    result<-"No nearby heterozygous SNPs to confirm"
  } else if(sum(phasing_info$n_ref_reads_against)>0 & sum(phasing_info$n_ref_reads_against)>(sum(phasing_info$n_ref_reads_supporting)*0.2)){
    result<-"Both alleles confirmed with reference allele"
  } else {
    result<-"May have biased allele sequencing or LOH - suggest further confirmation"
  }
  return(result)
}

create_mut_df=function(mut,tree,matrices) {
  #Create reference set of sample sets that form clades - used in assessing PVVs and MAVs
  all_clades=unique(tree$edge[,2])
  all_clade_samples=lapply(all_clades,function(node) getTips(node=node,tree=tree))
  
  mut_df<-dplyr::bind_rows(mapply(function(samples,clade) {return(data.frame(NV=sum(matrices$NV[mut,samples]),NR=sum(matrices$NR[mut,samples]),clades=clade))},samples=all_clade_samples,clade=all_clades,SIMPLIFY = FALSE))
  mut_df$pos_test<-apply(mut_df,1,function(x){(x[2]>=12 & (x[1]/x[2])>=0.25)|(x[2]>=8 & (x[1]/x[2])>=0.3)|(x[2]>=6 & (x[1]/x[2])>=0.4)})
  mut_df$neg_test<-apply(mut_df,1,function(x){x[1]==0 & (x[2])>=10})
  return(mut_df)
}

create_MAV_df=function(mut1,mut2,tree,matrices) {
  #Create reference set of sample sets that form clades - used in assessing PVVs and MAVs
  all_clades=unique(tree$edge[,2])
  all_clade_samples=lapply(all_clades,function(node) getTips(node=node,tree=tree))
  
  MAV_df<-dplyr::bind_rows(mapply(function(samples,clade) {return(data.frame(NV1=sum(matrices$NV[mut1,samples]),NV2=sum(matrices$NV[mut2,samples]),NR=(sum(matrices$NV[mut2,samples])+sum(matrices$NR[mut1,samples])),clades=clade))},samples=all_clade_samples,clade=all_clades,SIMPLIFY = FALSE))
  MAV_df$mut1_pos_test<-apply(MAV_df,1,function(x){(x[3]>=12 & (x[1]/x[3])>=0.2)|(x[3]>=7 & (x[1]/x[3])>=0.25)|(x[3]>=5 & (x[1]/x[3])>=0.5)})
  MAV_df$mut2_pos_test<-apply(MAV_df,1,function(x){(x[3]>=12 & (x[2]/x[3])>=0.2)|(x[3]>=7 & (x[2]/x[3])>=0.25)|(x[3]>=5 & (x[2]/x[3])>=0.5)})
  MAV_df$neg_test<-apply(MAV_df,1,function(x){sum(x[1:2])==0 & (x[3])>=10})
  
  return(MAV_df)
}

#Function to determine the "lesion path" and the fixed or "pure" subclades thrown off by the lesion
get_pure_subclades=function(mut1,mut2=NULL,lesion_node,tree,matrices) {
  if(is.null(mut2)) {test_type="PVV"} else {test_type="MAV"}
  print(paste("Testing",test_type))
  
  if(test_type=="PVV"){mut_df=create_mut_df(mut=mut1,tree=tree,matrices=matrices)} else {mut_df=create_MAV_df(mut1=mut1,mut2=mut2,tree=tree,matrices=matrices)}
  
  #1. get daughter nodes of lesion node
  lesion_children=get_node_children(lesion_node,tree=tree)
  if(length(lesion_children)>2) { #if initial_lesion_node is at site of polytomy, drop the negative branches of the polytomy
    print("Removing polytomy")
    if(test_type=="PVV") {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$pos_test[mut_df$clades%in%nodes]))})
    } else {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])|any(mut_df$mut2_pos_test[mut_df$clades%in%nodes]))})
    }
    lesion_children<-lesion_children[keep_children]
  }
  
  #Test these daughter nodes to see if they are "pure positive", "pure negative" or "mixed"
  if(test_type=="PVV") {types=get_node_types(lesion_children,mut_df,tree=tree)} else {types=get_MAV_node_types(lesion_children,mut_df,tree=tree)}
  if(sum(types=="mixed")>1) {
    stop(return("More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion"))
  } else if(length(unique(types))==1){
    stop(return(ifelse(test_type=="PVV","Not PVV","Not MAV")))  
  }
  
  names(lesion_children)<-types
  pure_subclades=lesion_children[names(lesion_children)!="mixed"]
  
  while(any(names(lesion_children)=="mixed")) {
    lesion_node=lesion_children["mixed"] #Get the new lesion node for this iteration (the "mixed" descendant of the previous lesion node)
    lesion_children=get_node_children(lesion_node,tree=tree)
    if(test_type=="PVV") {types=get_node_types(lesion_children,mut_df,tree=tree)} else {types=get_MAV_node_types(lesion_children,mut_df,tree=tree)}
    if(sum(types=="mixed")>1) {stop(return("More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion"))}
    names(lesion_children)=types
    if(!any(types=="mixed")) {
      pure_subclades=c(pure_subclades,lesion_children)
      break
    } else {
      lesion_node=lesion_children["mixed"]
      pure_subclades=c(pure_subclades,lesion_children[types!="mixed"])
    }
  }
  return(pure_subclades)
}

get_mixed_subclades=function(mut1,mut2=NULL,lesion_node,tree,matrices) {
  if(is.null(mut2)) {test_type="PVV"} else {test_type="MAV"}
  print(paste("Testing",test_type))
  
  if(test_type=="PVV"){mut_df=create_mut_df(mut=mut1,tree=tree,matrices=matrices)} else {mut_df=create_MAV_df(mut1=mut1,mut2=mut2,tree=tree,matrices=matrices)}
  
  #1. get daughter nodes of lesion node
  lesion_children=get_node_children(lesion_node,tree=tree)
  if(length(lesion_children)>2) { #if initial_lesion_node is at site of polytomy, drop the negative branches of the polytomy
    print("Removing polytomy")
    if(test_type=="PVV") {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$pos_test[mut_df$clades%in%nodes]))})
    } else {
      keep_children=sapply(lesion_children, function(node) {nodes=c(node,get_all_node_children(node,tree=tree)); return(any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])|any(mut_df$mut2_pos_test[mut_df$clades%in%nodes]))})
    }
    lesion_children<-lesion_children[keep_children]
  }
  
  #Test these daughter nodes to see if they are "pure positive", "pure negative" or "mixed"
  if(test_type=="PVV") {types=get_node_types(lesion_children,mut_df,tree=tree)} else {types=get_MAV_node_types(lesion_children,mut_df,tree=tree)}
  names(lesion_children)<-types
  if(sum(types=="mixed")>1) {
    stop(return("More than one mixed subclade identified - indicative that not caused by a persistent DNA lesion"))
  } else if(sum(types=="mixed")==0){
    stop(return("No mixed daughters"))
  } else {
    return(lesion_children["mixed"])
  }
}

get_file_paths_and_project=function(dataset,Sample_ID) {
  if(dataset=="MSC_fetal") {
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/Tree_",Sample_ID,".tree")
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/Filtered_mut_set_annotated_",Sample_ID)
    project=read.csv("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/MSC_BMT/Samples_project_reference.csv",header=T)
    project<-project[,c("Sample","Project")]
    colnames(project)<-c("sample","project")
    sex=NA
  } else if(dataset=="EM") {
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/tree_",Sample_ID,"_standard_rho01.tree")
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/annotated_mut_set_",Sample_ID,"_standard_rho01")
    project_ref=read.csv("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/EM/Samples_project_ref.csv",header=T)
    project_ref<-project_ref[,c(1,3)]
    colnames(project_ref)<-c("sample","project")
    sample=substr(Sample_ID,1,5)
    project=as.numeric(project_ref$project[project_ref$sample==sample])
    sex=NA
  } else if(dataset=="KY") {
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/",Sample_ID,"_rmix_consense_tree_no_branch_lengths_1811.tree") 
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/Filtered_muts_",Sample_ID)
    project_ref=read.csv("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/KY/Samples_project_ref_KY.csv",header=T)
    project=as.numeric(project_ref$project[project_ref$sample==Sample_ID])
    sex=NA
  } else if(dataset=="MSC_BMT") {
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/tree_",Sample_ID,"_m40_postMS_reduced_pval_post_mix.tree")
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/annotated_mut_set_",Sample_ID,"_m40_postMS_reduced_pval_post_mix")
    project=read.csv("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/MSC_BMT/Samples_project_reference.csv",header=T)
    project<-project[,c("Sample","Project")]
    colnames(project)<-c("sample","project")
    sex_vec=c(Pair11="male",Pair13="male",Pair21="male",Pair28="female",Pair31="male",Pair40="male")
    sex=sex_vec[Sample_ID]
  } else if(dataset=="PR") {
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/",Sample_ID,"/snp_tree_with_branch_length_polytomised.tree")
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/",Sample_ID,"/Filtered_muts_",Sample_ID)
    project=read.csv("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/PR/Samples_project_ref_PR.csv",header=T)
    project<-project[,c("sample","project")]
    sex=NA
  } else if(dataset=="MF"){
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/tree_",Sample_ID,"_noMixed.tree")
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/filtered_muts_",Sample_ID,"_noMixed")
    project=2305
    sex=NA
  } else if(dataset=="NW"){
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/tree_",Sample_ID,".tree")
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/filtered_muts_",Sample_ID)
    project=read.csv("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/NW/Samples_project_ref_NW.csv",header=T)
    project<-project[,c("Sample","Project")]
    colnames(project)<-c("sample","project")
    sex=NA
  } else if(dataset=="SN"){
    tree_file_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/tree_",Sample_ID,".tree")
    filtered_muts_path=paste0("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/",dataset,"/Filtered_muts_",Sample_ID)
    project=read.csv("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/SN/Samples_project_ref_SN.csv",header=T)
    project<-project[,c("sample","project")]
    colnames(project)<-c("sample","project")
    sex=NA
  }
  return(list(tree_file_path=tree_file_path,filtered_muts_path=filtered_muts_path,project=project,sex=sex))
}

#Combining the MAVs into a consensus Ref and Alt covering the same positions is complicated by the fact that they may report slightly different sequences (e.g. if one is an SNV and the other an MNV)
#This function is to define a reference set
establish_ref_and_alt=function(Ref1,Ref2,Alt1,Alt2,Pos1,Pos2) {
  if(Ref1==Ref2) {
    Ref<-Ref1 
  } else if(Pos1!=Pos2){
    Ref_vec1=as.vector(str_split(Ref1,"",simplify=T))
    names(Ref_vec1)=seq_along(Ref_vec1)+Pos1-1
    
    Ref_vec2=as.vector(str_split(Ref2,"",simplify=T))
    names(Ref_vec2)=seq_along(Ref_vec2)+Pos2-1
    
    Ref_vec<-c(Ref_vec1,Ref_vec2[!names(Ref_vec2)%in%names(Ref_vec1)])
    Ref<-paste0(Ref_vec,collapse="")
    
    #Do the same for the Alt1 vector
    Alt1_vec=as.vector(str_split(Alt1,"",simplify=T))
    names(Alt1_vec)=seq_along(Alt1_vec)+Pos1-1
    Alt1_new_vec=c(Alt1_vec,Ref_vec[!names(Ref_vec)%in%names(Ref_vec1)])
    Alt1<-paste0(Alt1_new_vec,collapse="")
    
    #Do the same for the Alt1 vector
    Alt2_vec=as.vector(str_split(Alt2,"",simplify=T))
    names(Alt2_vec)=seq_along(Alt2_vec)+Pos2-1
    Alt2_new_vec=c(Ref_vec[!names(Ref_vec)%in%names(Ref_vec2)],Alt2_vec)
    Alt2<-paste0(Alt2_new_vec,collapse="")
    
  } else if(nchar(Ref2)>nchar(Ref1)){
    Ref<-Ref2 #The longer ref is the "new ref"
    
    #Set this up as a vector named by the position of each base
    Ref_vec=as.vector(str_split(Ref,"",simplify=T))
    names(Ref_vec)=seq_along(Ref_vec)+Pos2-1
    
    #Do the same for the Alt1 vector
    Alt1_vec=as.vector(str_split(Alt1,"",simplify=T))
    names(Alt1_vec)=seq_along(Alt1_vec)+Pos1-1
    
    #Make the new Alt1 by replacing the matching positions of the Ref vec
    Alt1_new_vec<-Ref_vec
    Alt1_new_vec[names(Alt1_vec)]<-Alt1_vec
    Alt1<-paste0(Alt1_new_vec,collapse="")
    
  } else if(nchar(Ref1)>nchar(Ref2)){
    Ref<-Ref1 #The longer ref is the "new ref"
    
    #Set this up as a vector named by the position of each base
    Ref_vec=as.vector(str_split(Ref,"",simplify=T))
    names(Ref_vec)=seq_along(Ref_vec)+Pos1-1
    
    #Do the same for the Alt1 vector
    Alt2_vec=as.vector(str_split(Alt2,"",simplify=T))
    names(Alt2_vec)=seq_along(Alt2_vec)+Pos2-1
    
    #Make the new Alt1 by replacing the matching positions of the Ref vec
    Alt2_new_vec<-Ref_vec
    Alt2_new_vec[names(Alt2_vec)]<-Alt2_vec
    Alt2<-paste0(Alt2_new_vec,collapse="")
  }
  return(list(Ref=Ref,Alt1=Alt1,Alt2=Alt2))
}

#This looks through the output from the "Phase_MAVs.R" script and extracts a phasing summary
extract_MAV_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$positive_subclade_res)) {
    stop(return("No result"))
  } else {
    res<-list$positive_subclade_res
  }
  
  if(class(res)=="character") {
    stop(return(res))
  } else if(class(res)=="list"){
    res_vec=unlist(res)
  }
  
  if(length(res_vec)==1) {
    stop(return(res_vec))
  } else if(length(unique(res_vec))==1) {
    return(res_vec[1])
  } else {
    res_vec_MAV<-res_vec[grepl("pure_mut1",names(res_vec))&grepl("pure_mut2",names(res_vec))]
    if(length(res_vec_MAV)==0) {
      stop(return("Unable to confirm phasing"))
    } else if(any(res_vec_MAV=="Same phasing confirmed")) {
      return("Same phasing confirmed")
    } else {
      return("Unable to confirm phasing")
    }
  }
}

#A set of slightly fiddley functions to used in the functions to extract the phasing summary from the PVV phasing info
confirm_het_SNP=function(phasing_info_df){
  if(class(phasing_info_df)=="logical") {
    return(NULL)
  } else {
    true_het<-phasing_info_df$alt_phases_with_base!=phasing_info_df$ref_phases_with_base
    result=data.frame(SNP_site=phasing_info_df$SNP_site,res=sapply(true_het,function(res) {
      if(is.na(res)) {
        return(NA)
      } else if(res) {
        return("Heterozygous")
      } else {
        return("Not heterozygous")
      }
    }))
    return(result) 
  }
}

return_het_SNPs_from_positive_clades=function(positive_subclade_phasing_info) {
  res=lapply(positive_subclade_phasing_info,confirm_het_SNP)
  res[unlist(lapply(res,is.null))]<-NULL
  comb_res=Reduce(f=function(df1,df2) {return(full_join(df1,df2,by="SNP_site"))},res)
  if(!is.null(comb_res)&!all(is.na(comb_res[,grepl("res",colnames(comb_res))]))){
    final_res=apply(comb_res[,-1,drop=F],1,function(x) {res<-x[!is.na(x)]; if(length(unique(res))==1){return(res[1])}else{return(NA)}})
    het_SNPs=comb_res$SNP_site[final_res=="Heterozygous" & !is.na(final_res)]
    if(length(het_SNPs)>0) {
      return(het_SNPs) 
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

get_alt_base=function(SNP,positive_subclade_phasing_info) {
  alts=unlist(lapply(positive_subclade_phasing_info,function(df) {
    if(class(df)=="logical") {
      return(NA)
    } else {
      return(df$alt_phases_with_base[df$SNP_site==SNP])
    }
  }))
  alts<-alts[!is.na(alts)]
  if(length(unique(alts))>1) {
    return("Conflicting results")
  } else {
    names(alts)<-NULL
    return(alts[1])
  }
}

assess_presence_of_alt_allele=function(alt_bases,negative_subclade_phasing_info) {
  het_SNP_sites=names(alt_bases)
  res=lapply(negative_subclade_phasing_info,function(df) {
    if(class(df)=="logical") {
      return("Alt allele not confirmed")
    } else if(!any(het_SNP_sites%in%df$SNP_site)) {
      return("Alt allele not confirmed")
    } else {
      res2=sapply(het_SNP_sites,function(SNP) {
        alt_base=alt_bases[SNP]
        if(SNP%in%df$SNP_site){
          if(df$ref_phases_with_base[df$SNP_site==SNP]==alt_base) {
            return("Alt allele reads present")
          } else if(df$ref_phases_with_base[df$SNP_site==SNP]!=alt_base & df$n_ref_reads_against[df$SNP_site==SNP]>0) {
            return("Alt allele reads present")
          } else {
            return("Alt allele not confirmed")
          } 
        } else {
          return(NA)
        }
      })
      res2<-res2[!is.na(res2)]
      if(any(res2=="Alt allele reads present" )) {
        return("Alt allele reads present")
      } else {
        max_reads=max(df$n_ref_reads_supporting[df$SNP_site%in%het_SNP_sites])
        return(paste("Alt allele not confirmed with maximum of",max_reads,"reads supporting the other allele"))
      } 
    }
  })
  return(res)
}

#This function extracts the positive subclade phasing info results from the results list
extract_MAV_pos_clade_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$positive_subclade_res)) {
    stop(return("No result"))
  } else {
    res<-list$positive_subclade_res
  }
  
  if(class(res)=="character") {
    stop(return(res))
  } else if(class(res)=="list"){
    res_vec=unlist(res)
  }
  
  if(length(res_vec)==1) {
    stop(return(res_vec))
  } else if(length(unique(res_vec))==1) {
    return(res_vec[1])
  } else {
    res_vec_MAV<-res_vec[grepl("pure_mut1",names(res_vec))&grepl("pure_mut2",names(res_vec))]
    if(length(res_vec_MAV)==0) {
      stop(return("Unable to confirm phasing"))
    } else if(any(res_vec_MAV=="Same phasing confirmed")) {
      return("Same phasing confirmed")
    } else {
      return("Unable to confirm phasing")
    }
  }
}

#This function extracts the positive subclade phasing info results from the results list
extract_PVV_pos_clade_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$positive_subclade_res)) {
    stop(return("No result"))
  } else {
    res_pos<-list$positive_subclade_res
  }
  
  if(class(res_pos)=="character") {
    stop(return(res_pos))
  } else if(class(res_pos)=="list"){
    res_pos_vec=unlist(res_pos)
  }
  
  if(length(res_pos_vec)==1) {
    stop(return(res_pos_vec))
  } else if(length(unique(res_pos_vec))==1) {
    return(res_pos_vec[1])
  } else {
    if(any(res_pos_vec=="Same phasing confirmed")) {
      return("Same phasing confirmed in at least one subclade")
    } else if(any(res_pos_vec=="Non-matching phasing confirmed")) {
      return("Non-matching phasing confirmed in at least one subclade")
    } else {
      return("Unable to confirm phasing")
    }
  }
}

#This function examines the read counts of the negative subclades of the PVV to see if they include reads that match
#the phasing of the mutant allele in the positive subclades. If they do this means (1) there is no LOH, (2) both alleles have been sequenced
extract_PVV_neg_clade_phasing_summary=function(list) {
  if(class(list)!="list") {
    stop(return("No result"))
  } else if(is.null(list$negative_subclade_res)) {
    stop(return("No result"))
  } else {
    res_neg<-list$negative_subclade_res
  }
  
  if(class(res_neg)=="character") {
    stop(return(res_neg))
  } else if(class(res_neg)=="list"){
    res_neg_vec=unlist(res_neg)
  } else if(is.na(res_neg)) {
    stop(return(NA))
  }
  
  if(length(res_neg_vec)==1) {
    res<-res_neg_vec
  } else if(length(unique(res_neg_vec))==1) {
    res<-res_neg_vec[1]
  } else if(any(res_neg_vec=="Both alleles confirmed with reference allele")){
    res<-"Both alleles confirmed with reference allele in at least one subclade"
  } else {
    res<-res_neg_vec
  } 
  
  if(any(res=="May have biased allele sequencing or LOH - suggest further confirmation")) {
    pos_clades=which(names(list$phasing_info_by_subclade)%in%c("pure_positive","pure_mut1","pure_mut2"))
    het_SNPs=return_het_SNPs_from_positive_clades(list$phasing_info_by_subclade[pos_clades])
    #print(het_SNPs)
    if(!is.null(het_SNPs)) {
      alt_bases<-sapply(het_SNPs,function(SNP) {get_alt_base(SNP,list$phasing_info_by_subclade[pos_clades])})
      if(all(alt_bases=="Conflicting results")) {
        stop(return("Positive clades have non-matching phasing"))
      }
      alt_bases<-alt_bases[!alt_bases=="Conflicting results"]
      neg_clades=which(names(list$phasing_info_by_subclade)=="pure_negative")
      res<-unlist(assess_presence_of_alt_allele(alt_bases,list$phasing_info_by_subclade[neg_clades]))
      res<-unique(res)
      if(length(res)>1) {
        if(any(res=="Alt allele reads present")) {
          res<-"Alt allele reads present in at least one subclade"
        } else {
          res<-paste(res,collapse=",")
        }
      }
    } else {
      res<-"No nearby heterozygous SNPs to confirm"
    }
  }
  return(res) 
}

#ASCAT copy number functions
get_cn=function(cn_summary_file){
  ##cat("opening ",cn_summary_file,"\n")
  if(!file.exists(cn_summary_file)){
    warning(sprintf("%s: does not exist",cn_summary_file))
    return(NULL)
  }
  cn=read.csv(cn_summary_file,header = FALSE)
  cn$start=cn$V3
  cn$end=cn$V4
  cn$chr=cn$V2
  ###V7=total copy number
  ## V8=minor allele copy number
  cn$major=cn$V7-cn$V8
  cn$minor=cn$V8
  cn[,-grep("^V",colnames(cn))]
}
get_ASCAT_minor_allele_cn=function(Chrom,Pos,sample,project){
  if(is.data.frame(project)) {
    sample_project<-project$project[project$sample==sample]<-project$project[project$sample==sample]
  } else {
    sample_project<-project
  }
  file = paste0("/nfs/cancer_ref01/nst_links/live/", sample_project, "/", sample, "/", sample, ".ascat_ngs.summary.csv")
  cn=get_cn(file)
  if(!is.null(cn)) {
    minor_allele_cn=cn$minor[cn$chr==Chrom & cn$start<Pos & cn$end>Pos]
    return(minor_allele_cn)
  } else {
    return(NA)
  }
}

get_mean_ASCAT_minor_allele_cn=function(Chrom,Pos,samples,project) {
  cn_vec=sapply(samples, function(sample) {
    #print(sample)
    cn=get_ASCAT_minor_allele_cn(Chrom = Chrom,Pos=Pos,sample=sample,project=project)
    return(cn)
  })
  #print(cn_vec)
  return(mean(cn_vec,na.rm = T))
}

##FUNCTIONS FOR THE ANALYSIS OF LESION SEGREAGATION DATA
#estimate the parameters pf
estimate_gamma_params=function(value_vec,log_rate_range=c(-2,1),shape_range=c(1,5)) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rate_vec = 10^(seq(log_rate_range[1],log_rate_range[2],by=0.1)) # rho will be bounded within 1e-6 and 0.89
  shape_vec=seq(shape_range[1],shape_range[2],0.05)
  params_grid=expand_grid(rate_vec,shape_vec)
  ll = sapply(1:nrow(params_grid), function(i) {shape=params_grid$shape_vec[i]; rate=params_grid$rate_vec[i];sum(dgamma(x=value_vec, shape=shape,rate=rate,log = T))})
  return(params_grid[which.max(ll),])
}

#Updated version of the squash_tree function that allows you to squash from the root, as well as from the tips
squash_tree=function(tree,cut_off=50,from_root=F) {
  if(from_root){
    idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
    new_edge_lengths=nodeHeights(tree)[idxs_to_squash,2]-cut_off #work-out the edge lengths that these should be such that they finish at the cut-off
    tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
    
    tree$edge.length[nodeHeights(tree)[,2]<=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
    return(tree)
  } else {
    tree$edge.length[nodeHeights(tree)[,1]>=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
    idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
    new_edge_lengths=cut_off - nodeHeights(tree)[idxs_to_squash,1] #work-out the edge lengths that these should be such that they finish at the cut-off
    tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
    return(tree)
  }
}

#First version of the sharedness stat - as per NW. However, this will give higher values of sharedness with smaller trees
calculate_sharedness_stat=function(tree) {
  prop_samples<-sapply(tree$edge[,2],function(node) {
    prop_samples<-length(getTips(tree,node))/length(tree$tip.label)
    return(prop_samples)
  })
  mean_w<-weighted.mean(x=prop_samples,w=tree$edge.length)
  return(mean_w)
}

#Second version of the sharedness stat. Minus 1 from the numerator & denominator.
calculate_sharedness_stat_2=function(tree) {
  prop_samples<-sapply(tree$edge[,2],function(node) {
    prop_samples<-(length(getTips(tree,node))-1)/(length(tree$tip.label)-1)
    return(prop_samples)
  })
  mean_w<-weighted.mean(x=prop_samples,w=tree$edge.length)
  return(mean_w)
}

#Count the number of internal nodes above a certain height i.e. for calculating the "post-developmental nodes"
count_internal_nodes=function(tree,cut_off=50){
  nodeheights=nodeHeights(tree)
  internal_nodes_above_cutoff=sum(nodeheights[,2]>cut_off & !tree$edge[,2]%in%1:length(tree$tip.label))
  return(internal_nodes_above_cutoff)
}