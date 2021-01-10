#This function is required in the filtering function
get_ancestor_node=function(node,tree,degree=1){ #to get the 1st degree ancestor (i.e. the direct parent) use degree=1.  Use higher degrees to go back several generations.
  curr<-node
  for(i in 1:degree){
    curr=tree$edge[which(tree$edge[,2]==curr),1]
    if(curr==(1+length(tree$tip.label))) {stop(return(curr))}
  }
  return(curr)
}

#Define the "get_node_types" function required for following the lesion journey in the case of PVVs
#It uses the mutation dataframe ("mut_df") to work out whether daughter branches of a node are (a) a mutant allele (b) wild-type or (c) mixed
get_node_types=function(lesion_children,mut_df) {
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

#As above, but can incorporate two alternative mutant alleles, therefore suitable for MAVs
get_MAV_node_types=function(lesion_children,mut_df) {
  types=sapply(lesion_children, function(node) {
    nodes=c(node,get_all_node_children(node,tree))
    if(all(!mut_df$neg_test[mut_df$clades%in%nodes])&any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])&!any(mut_df$mut2_pos_test[mut_df$clades%in%nodes])){
      return("pure_mut1")
    } else if(all(!mut_df$neg_test[mut_df$clades%in%nodes])&any(mut_df$mut2_pos_test[mut_df$clades%in%nodes])&!any(mut_df$mut1_pos_test[mut_df$clades%in%nodes])){
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
      near_muts=which(details_by_chrom$Pos>Pos&
                        details_by_chrom$Pos<=(Pos+region_size)&
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
    new_df=Reduce(rbind,new_df)
    
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
get_multi_allelic_variant_list=function(details) {
  Chroms=c(1:22,"X","Y")
  details$Ref=as.character(details$Ref)
  details$Alt=as.character(details$Alt)
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
}


