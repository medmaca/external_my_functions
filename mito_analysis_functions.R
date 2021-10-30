#Plotting function for this script - overlays the VAFs of mitochondrial mutations across the tree (colour scale)
#Rescales them to be between the minimum & maximum VAF of the samples (for maximum contrast)
plot_VAF=function(tree,
                  details,
                  matrices,
                  node,
                  mut1,
                  colours=c("dark gray","red"),
                  min_vaf=0,
                  max_vaf=1,
                  cex=0.4,
                  #show_pval=FALSE,
                  ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  
  mut_colfunc = colorRampPalette(colours)
  mut_colscale = mut_colfunc(101)
  
  #Get the vaf
  samples=getTips(tree = tree,node=node)
  mut1_reads=sum(matrices$NV[mut1,samples])
  total_reads=sum(matrices$NR[mut1,samples])
  
  VAF=mut1_reads/total_reads
  
  newrange=c(0,1)
  xrange<-range(min_vaf,max_vaf)
  mfac <- (newrange[2] - newrange[1])/(xrange[2] - xrange[1])
  plot_x=newrange[1] + (VAF - xrange[1]) * mfac
  
  branch_col=mut_colscale[1+100*round(plot_x,digits=2)] #Now how "concentrated" should the mut colour be
  
  info=get_edge_info(tree,details,node)
  
  #Plot the branches using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=branch_col,lend=1,...)
  }
  #Print the mutation name (only print it once)
  if(node==1) {text(1,1,pos=4,paste(mut1,": Maximum VAF is",round(max_vaf,digits = 3),"Minimum VAF is:",round(min_vaf,digits = 3)))}
}

plot_multi_VAF=function(tree,
                        details,
                        matrices,
                        node,
                        muts,
                        colours=c("#4081ec", "#17324d", "#65b5d8", "#21247b", "#92a654", "#02531d", "#68c030", "#5e132a", "#ec929b", "#965d39", "#0fb381", "#fd2c3b", "#fe8f06", "#b70d61", "#087394"),
                        cex=0.4,
                        ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  
  if(length(muts)>length(colours)) {
    mut_cols<-colorRampPalette(colours)(length(muts))
  } else {
    mut_cols<-colours[1:length(muts)]
  }
  names(mut_cols)<-muts
  
  #Get the vaf
  samples=getTips(tree = tree,node=node)
  mut_reads=rowSums(matrices$NV[muts,samples,drop=F])
  total_reads=rowSums(matrices$NR[muts,samples,drop=F])
  
  VAF=mut_reads/total_reads
  VAF<-VAF[VAF>0.05]
  if(length(VAF)==0) {
    branch_col<-"lightgrey"
  } else if(length(VAF)==1) {
    branch_col<-mut_cols[names(VAF)]
  } else if(length(VAF)>1) {
    top2=names(sort(VAF,decreasing=T))[1:2]
    cols=mut_cols[top2]
    
    mut_colfunc = colorRampPalette(cols)
    mut_colscale = mut_colfunc(101)
    branch_col=mut_colscale[1+round(100*VAF[top2][1]/sum(VAF[top2]))] #Now how "concentrated" should the mut colour be
  }
  
  info=get_edge_info(tree,details,node)
  
  #Plot the branches using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=branch_col,lend=1,...)
  }
  #Print the mutation name (only print it once)
  #if(node==1) {text(1,1,pos=4,paste(mut1,": Maximum VAF is",round(max_vaf,digits = 3),"Minimum VAF is:",round(min_vaf,digits = 3)))}
}

reverse_germline=function(matrices) {
  mean_vaf=rowMeans(matrices$NV/matrices$NR)
  matrices$NV[mean_vaf>0.5,]<-(matrices$NR[mean_vaf>0.5,]-matrices$NV[mean_vaf>0.5,])
  matrices$SW[mean_vaf>0.5,]<-0 #Shearwater not set up to call somatic reversion mutations, therefore just set these to 0 across samples
  return(matrices)
}

generate_mito_matrices=function(PD_number,tree_file_path,pileup_folder=NULL,shearwater_calls_file=NULL,haplocheck_calls_folder=NULL,reverse_germline=T,run_bb=F){
  library(dplyr)
  library(stringr)
  tree=read.tree(tree_file_path)
  
  #Read in the pileup for a given pair
  if(!is.null(pileup_folder)){
    files=list.files(pileup_folder,pattern = PD_number,full.names = T)
    if(length(files)==0) {print("No pile up files found for this PD_number in the specified pileup folder");stop(return(NULL))}
    temp=read.csv(files[1])
    mt_dat=array(0,dim=c(nrow(temp),ncol(temp),length(files)))
    genotypes=c("A","T","C","G","DEL","INS","a","t","c","g","del","ins")
    sample_names=gsub(pattern = "_MT_count.csv",replacement = "", x = list.files(pileup_folder,pattern = PD_number,full.names = F))
    dimnames(mt_dat)=list(1:nrow(temp),genotypes,sample_names)
    cat("Reading in pile up files\n")
    for(i in 1:length(files)) {
      file=files[i]
      import=read.csv(file)
      colnames(import)=genotypes
      mt_dat[,,i]<-as.matrix(import)
    }
  }
  
  #Read in the shearwater calls for this individual
  if(!is.null(shearwater_calls_file)){
    cat("Reading in shearwater calls file\n")
    shearwater_calls=read.delim(shearwater_calls_file)
    shearwater_calls$pos=as.character(shearwater_calls$pos)
    shearwater_calls<-shearwater_calls[shearwater_calls$pos!="3107",]
    shearwater_calls$mut_ref=apply(shearwater_calls[,c("chr","pos","ref","mut")],1,paste,collapse="_")
    shearwater_calls$mut_ref<-str_replace(shearwater_calls$mut_ref,pattern="-",replacement="DEL")
    shearwater_calls<-shearwater_calls%>%filter(sampleID%in%tree$tip.label)
  } else {
    shearwater_calls<-NULL
  }
  
  #Read in the haplocheck calls for this individual
  if(!is.null(haplocheck_calls_folder)){
    cat("Checking for matching files in the haplocheck calls folder\n")
    files=list.files(haplocheck_calls_folder,pattern = PD_number,full.names = T)
    if(length(files)==0) {
      cat("No haplocheck files found for this PD_number in the specified haplocheck folder\n");stop(return(NULL))
    } else {
      cat(paste(length(files),"matching files found in the haplocheck calls folder\n"))
      }
    temp=read.delim(files[1],stringsAsFactors = F)
    sample_names=gsub(pattern = ".txt",replacement = "", x = list.files(haplocheck_calls_folder,pattern = PD_number,full.names = F))
    cat("Reading in haplocheck calls file\n")
    haplocheck_calls<-dplyr::bind_rows(Map(Sample=sample_names,file=files,function(Sample,file){
      sample_haplo<-read.delim(file,stringsAsFactors = F)
      sample_haplo$sampleID<-sampleID
      return(sample_haplo)
    }))
  } else {
    haplocheck_calls<-NULL
  }
  
  
  #Only include tips on the tree that have had mitochondrial mutation calls done.
  tree<-drop.tip(tree,tip = tree$tip.label[!tree$tip.label%in%sample_names])
  
  #Get each mutation call that has passed Shearwater in at least 1 sample
  mut_table=table(shearwater_calls$mut_ref)
  mut_table_shared=mut_table[mut_table>1]#Select those that are positive in >2 samples
  mut_all=sort(unique(shearwater_calls$mut_ref))
  
  #Create NV and NR matrices.
  cat("Creating the NV and NR matrices\n")
  NV=NR=matrix(0,nrow=length(mut_all),ncol=length(tree$tip.label))
  dimnames(NV)=dimnames(NR)=list(mut_all,tree$tip.label)
  for(i in 1:length(mut_all)) {
    if(i%%1000==0) {print(i)}
    mut=mut_all[i]
    pos=as.numeric(str_split(mut,pattern="_",simplify=T)[,2])
    ref=str_split(mut,pattern="_",simplify=T)[,3]
    alt=str_split(mut,pattern="_",simplify=T)[,4]
    dep_counts<-apply(mt_dat[as.character(pos),,],2,sum)
    mut_counts<-apply(mt_dat[as.character(pos),c(alt,tolower(alt)),],2,sum)
    NR[i,]<-dep_counts[tree$tip.label]
    NV[i,]<-mut_counts[tree$tip.label]
  }
  
  #Convert to df, keeping only those columns that are samples in the tree
  NV=as.data.frame(NV[,tree$tip.label])
  NR=as.data.frame(NR[,tree$tip.label])
  
  SW_mat<-matrix(0,nrow=nrow(NV),ncol=ncol(NV),dimnames=dimnames(NV))
  for(i in 1:nrow(shearwater_calls)){
    SW_mat[shearwater_calls$mut_ref[i],shearwater_calls$sampleID[i]]<-1
    }
  
  #Calculate over-dispersion. Most likely to be informative lineage markers.
  if(run_bb){
    res=beta.binom.filter(COMB_mats = list(NV=NV,NR=NR))
  } else {
    res<-NA
  }
  
  
  #Now add a global (aggregated counts for all samples) column
  NV$global=rowSums(NV)
  NR$global=rowSums(NR)
  
  matrices=list(NV=NV,NR=NR,SW=SW_mat)
  
  #If 'reverse germline' option selected, reverse the mut/ wt calls for those mutations that are more common than the wild type
  #(i.e. likely to have been mutant in the oocyte)
  if(reverse_germline){
    matrices=reverse_germline(matrices)
  }
  
  vaf=calculate_vaf(matrices$NV,matrices$NR)
  return(list(matrices=list(vaf=vaf,NV=NV,NR=NR,SW=SW_mat),rho_vals=res,tree=tree,sample_shearwater_calls=shearwater_calls,sample_haplocheck_calls=haplocheck_calls))
}

#Tree to find the latest (smallest) clade enclosing all the positive samples
find_latest_acquisition_node=function(tree,pos_samples) {
  curr_node<-which(tree$tip.label==pos_samples[1])
  while(!all(pos_samples%in%getTips(tree,curr_node))){
    curr_node<-getAncestors(tree,curr_node,type="parent")
  }
  return(curr_node)
}

#Adjust the add_heatmap function to allow large heatmap under the tree
add_mito_mut_heatmap=function(tree,heatmap,heatvals=NULL,border="white",heatmap_bar_height=0.05,cex.label=2){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.01*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-heatmap_bar_height*ymax
    #bot=top-(0.05/dim(heatmap)[1])*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border,lwd = 0.25)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=-0.5,y=0.5*(top+bot),pos = 2,cex = cex.label)
    }
    top=bot
  }
  tree
}


#####################################################################################
# FUNCTION
#####################################################################################

trinucleotide_plot = function (mutations, file_name=NULL, analysis_type, analysis_region) {
  list.of.packages <- c("BiocManager", "reshape2", "stringr", "readr", "tidyverse", "ggpubr", "SummarizedExperiment", "Rsamtools")
  suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
  
  # subset the mutations to the respective columns
  mutations <- unique(mutations[,c("chr","pos","ref","mut","donor")])
  mutations <- mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")),]
  
  if(analysis_region == "coding"){
    mutations <- mutations[which(mutations$pos %in% coding_region),]
  }else if(analysis_region == "d_loop"){
    mutations <- mutations[which(mutations$pos %in% d_loop_region),]
  }else if(analysis_region != "all_mtDNA"){
    mutations <- NULL
  }
  
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))
  
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs_heavy = table(paste(mutations$sub[which(mutations$ref %in% c("A","G"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],3,3),sep="-"),sep=","))
  freqs_light = table(paste(mutations$sub[which(mutations$ref %in% c("C","T"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],3,3),sep="-"),sep=","))
  
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_heavy_full = freqs_heavy[full_vec]; freqs_heavy_full[is.na(freqs_heavy_full)] = 0; names(freqs_heavy_full) = full_vec
  freqs_light_full = freqs_light[full_vec]; freqs_light_full[is.na(freqs_light_full)] = 0; names(freqs_light_full) = full_vec
  
  
  if(analysis_type == "obs_exp"){
    if(analysis_region == "coding"){
      heavy_base_freqs <- (mtdna_trinuc_freq["coding_heavy",] / sum(mtdna_trinuc_freq["coding_heavy",]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- (mtdna_trinuc_freq["coding_light",] / sum(mtdna_trinuc_freq["coding_light",]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
      
    }else if(analysis_region == "d_loop"){
      heavy_base_freqs <- (mtdna_trinuc_freq["d_loop_heavy",] / sum(mtdna_trinuc_freq["d_loop_heavy",]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- (mtdna_trinuc_freq["d_loop_light",] / sum(mtdna_trinuc_freq["d_loop_light",]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
      
    }else if(analysis_region == "all_mtDNA"){
      heavy_base_freqs <- colSums(mtdna_trinuc_freq[c("coding_heavy","d_loop_heavy"),]) / sum(colSums(mtdna_trinuc_freq[c("coding_heavy","d_loop_heavy"),]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- colSums(mtdna_trinuc_freq[c("coding_light","d_loop_light"),]) / sum(colSums(mtdna_trinuc_freq[c("coding_light","d_loop_light"),]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
      
    }else{
      freqs_heavy_full = NULL
      freqs_light_full = NULL
    }
  }
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  #dev.new(width=10,height=4)
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y_heavy = freqs_heavy_full; y_light = freqs_light_full; maxy = max(c(y_heavy,y_light))
  
  if(analysis_type == "obs_exp"){
    ylab = "Mutation frequency (Obs/Exp)"
  }else{
    ylab = "Mutation count"
  }
  
  h_heavy = barplot(y_heavy, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab)
  h_light = barplot(-y_light, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab, add = T)
  
  segments(y0 = maxy*1.5, y1 = maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
  segments(y0 = -maxy*1.5, y1 = -maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
  segments(y0 = 0, y1 = 0, x0 = 0.5, x1 = 192.5,  col = "black")
  abline(v = 0.5, col = "black")
  abline(v = 32.5, col = "black")
  abline(v = 64.5, col = "black")
  abline(v = 96.5, col = "black")
  abline(v = 128.5, col = "black")
  abline(v = 160.5, col = "black")
  abline(v = 192.5, col = "black")
  
  
  for (j in 1:length(sub_vec)) {
    xpos = h_heavy[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.25, xpos[2]+0.5, maxy*1.15, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.15, pos=3, labels=sub_vec[j])
  }
  if(!is.null(file_name)){
    dev.copy(pdf,file_name,width=12,height=5)
    dev.off()
  }
  #dev.off()
}

