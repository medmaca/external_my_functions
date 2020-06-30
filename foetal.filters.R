#Function to subset the list of matrices
list_subset = function(list, select_vector) {
  for (i in 1:length(list)) {
    list[[i]] = list[[i]][select_vector,]
  }
  return(list)
}

#Function to create pval matrix based on the NV and NR matrix
pval_matrix = function(mat, NV, NR, gender) {
  pval_mat <- matrix(0, nrow = nrow(NV), ncol = ncol(NV))
  if(gender == "male") {
    for(i in 1:nrow(NV)) {
    for (j in 1:ncol(NR)) {
      if (!mat$Chrom[i] %in% c("X","Y")) {pval_mat[i,j] <- binom.test(NV[i,j], NR[i,j], p = 0.5, alternative = "less")$p.value}
      else {pval_mat[i,j] <- binom.test(NV[i,j], NR[i,j], p = 0.95, alternative = "less")$p.value}
    }
    if (i %% 1000 == 0) {print(i)}
    }
  } else if(gender == "female") {
    for(i in 1:nrow(NV)) {
      for (j in 1:ncol(NR)) {
        pval_mat[i,j] <- binom.test(NV[i,j], NR[i,j], p = 0.5, alternative = "less")$p.value
      }
      if (i %% 1000 == 0) {print(i)}
    }
  }
  return(pval_mat)
}

#mat object needs the Chrom column with chromosome.  Returns the pval vector for each mutation.
germline.binomial.filter = function(mat, NV, NR, gender){
  XY_chromosomal = mat$Chrom%in%c("X","Y")
  autosomal = !XY_chromosomal
  
  if(gender=="female"){
    NV_vec = rowSums(NV) #create vector of COMBINED variant reads across all samples
    NR_vec = rowSums(NR) #create vector of COMBINED depth across all samples
    pval = rep(1,length(NV_vec))
    #For loop to calculate whether each one is likely to have come from a binomial distribution where true probability is 0.5 (i.e. would be expected for autosomal heterozygous germline variant)
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      } 
      if (n%%1000==0){
        print(n)
      }
    }
  }
  # IF MALE - need to separate off the XY chromosomes, as expected probability of germline variant is now close to 1.
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    NV_vec = rowSums(NV)[autosomal]
    NR_vec = rowSums(NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    NV_vec = rowSums(NV)[XY_chromosomal]
    NR_vec = rowSums(NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  return(pval)
}

#THE BETA-BINOMIAL FUNCTIONS
require(VGAM)
estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV,cutoff=0.3, binom.pval=F,pval.cutoff=0.05){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec=as.numeric(NR[k,]))
    if(binom.pval){
      mu = sum(as.numeric(NV[k,]))/sum(as.numeric(NR[k,]))
      LL0 = sum(dbinom(as.numeric(NV[k,]),as.numeric(NR[k,]),prob=mu))
      LL1 = sum(dbetabinom(as.numeric(NV[k,]),as.numeric(NR[k,]),prob=mu,rho=rho_est[k]))
      pval[k] = (1-pchisq(2*(LL1-LL0),df=1)) / 2
    }
    if (k%%1000==0){
      print(k)
    }
  }
  if(binom.pval){
    qval=p.adjust(pval,method="BH")
    flt_rho=qval<=pval.cutoff&rho_est>cutoff
  }else{
    flt_rho=rho_est>=cutoff
  }
  return(rho_est)
}

#FILTER VARIANTS WITH LOW VAF AMONGST POSITIVE SAMPLES
#Filters mutations present in multiple samples that are (on aggregate) unlikely to have come from true somatic mutation distribution (i.e. 0.5 in auto, 1 in XY)
low_vaf_in_pos_samples = function(mat, NV, NR, define_pos = 2) {
  if(all(c(nrow(mat) == nrow(NV),dim(NV) == dim(NR)))) {print("Input matrices are of correct dimensions")}
  pval=rep(0,nrow(mat))
  for(k in 1:nrow(mat)) {
    NV_vec <- NV[k,]
    NR_vec <- NR[k,]
    if(any(NV_vec >= define_pos)){
      NV_vec_pos <- NV_vec[which(NV_vec >= define_pos)]
      NR_vec_pos <- NR_vec[which(NV_vec >= define_pos)]
      if (mat$Chrom[k] %in% c("X","Y")) {
        pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
      }
      else {
        pval[k] <- binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
      }
      if(k %% 1000 == 0) {print(k)}
    }
  }
  return(pval)
}

#Removing columns with low coverage, or that are otherwise unwanted
remove_low_coverage_samples = function(COMB_mats,
                                       filter_params,
                                       min_sample_mean_cov,
                                       other_samples_to_remove = NULL,
                                       min_variant_reads_auto = 3,
                                       min_variant_reads_xy = 2) {
  mean_sample_cov = colMeans(COMB_mats$NR)
  low_cov_names = gsub(pattern = "_DEP", replacement = "", x = names(mean_sample_cov)[mean_sample_cov < min_sample_mean_cov])
  remove_cols <- which(gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats$NV)) %in% c(low_cov_names, other_samples_to_remove))
  cat("Removing",length(remove_cols),"samples\n")
  COMB_mats$NV <- COMB_mats$NV[,-remove_cols]
  COMB_mats$NR <- COMB_mats$NR[,-remove_cols]
  COMB_mats$PVal <- COMB_mats$PVal[,-remove_cols]
  
  null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto|(COMB_mats$NV >= min_variant_reads_xy & COMB_mats$mat$Chrom %in% c("X","Y"))) == 0
  cat(sum(null_remove),"mutations removed as no positives in any remaining samples.\n")
  COMB_mats = list_subset(COMB_mats, select_vector = !null_remove)
  filter_params = filter_params[!null_remove,]
  output = list(COMB_mats=COMB_mats,filter_params=filter_params)
  return(output)
}

#Functions for filtering from the filter_params and COMB_mats object, setting the desired cut-offs
assess_mean_depth = function(i, mat_list, AUTO_low_depth_cutoff, AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff) {
  if(!mat_list$mat$Chrom[i] %in% c("X","Y")) {
    ifelse((filter_params$mean_depth[i] >= AUTO_low_depth_cutoff & filter_params$mean_depth[i] <= AUTO_high_depth_cutoff), 1,0)
  } else {
    ifelse(filter_params$mean_depth[i] >= XY_low_depth_cutoff & filter_params$mean_depth[i] <= XY_high_depth_cutoff, 1,0)
  }
}

assess_max_depth_in_pos = function(i, mat_list, min_depth_auto, min_depth_xy) {
  if(!mat_list$mat$Chrom[i] %in% c("X","Y")) {
    ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_auto, 1,0)
  } else {
    ifelse(filter_params$max_depth_in_pos_samples[i] >= min_depth_xy, 1,0)
  }
}

get_filtered_mut_set = function(input_set_ID,
								COMB_mats,
                                filter_params,
                                gender,
                                retain_muts = NULL,
                                germline_pval_cutoff = -10,
                                rho_cutoff = 0.15,
                                XY_low_depth_cutoff = NULL,
                                XY_high_depth_cutoff = NULL,
                                AUTO_low_depth_cutoff=NULL,
                                AUTO_high_depth_cutoff=NULL,
                                pval_cutoff_dp2=NULL,
                                pval_cutoff_dp3=0.01,
                                min_depth_auto = 6,
                                min_depth_xy = 4,
                                min_pval_for_true_somatic = 0.1,
                                min_variant_reads_SHARED = 2,
                                min_pval_for_true_somatic_SHARED = 0.05) {
  #APPLY THE FILTERS - SET-UP THE "filter_pass" DATAFRAME
  filter_pass = filter_params
  filter_pass$germline_pval = ifelse(log10(filter_params$germline) < germline_pval_cutoff,1, 0)
  filter_pass$bb_rhoval = ifelse(filter_params$bb_rhoval > rho_cutoff, 1, 0)
  if(any(is.null(c(XY_low_depth_cutoff,XY_high_depth_cutoff,AUTO_low_depth_cutoff,AUTO_high_depth_cutoff)))) {
    sensitivity = 3 #number of standard deviations from mean depth to keep
    cat("Assigning mean-depth cutoffs automatically: removing mutations with a mean depth of >",sensitivity,"standard deviations from the mean\n")
    if(gender == "male") {
      mean_xy = mean(rowMeans(COMB_mats$NR)[COMB_mats$mat$Chrom %in% c("X","Y")])
      sd_xy = sd(rowMeans(COMB_mats$NR)[COMB_mats$mat$Chrom %in% c("X","Y")])
      if(is.null(XY_low_depth_cutoff)){XY_low_depth_cutoff = mean_xy - sensitivity*sd_xy}
      if(is.null(XY_high_depth_cutoff)) {XY_high_depth_cutoff = mean_xy + sensitivity*sd_xy}
      
      mean_auto = mean(rowMeans(COMB_mats$NR)[!COMB_mats$mat$Chrom %in% c("X","Y")])
      sd_auto = sd(rowMeans(COMB_mats$NR)[!COMB_mats$mat$Chrom %in% c("X","Y")])
      if(is.null(AUTO_low_depth_cutoff)) {AUTO_low_depth_cutoff = mean_auto - sensitivity*sd_auto} 
      if(is.null(AUTO_high_depth_cutoff)) {AUTO_high_depth_cutoff  = mean_auto +sensitivity*sd_auto}
      
    } else if(gender == "female") {
      mean_auto = mean(rowMeans(COMB_mats$NR))
      sd_auto = sd(rowMeans(COMB_mats$NR))
      if(is.null(AUTO_low_depth_cutoff)) {AUTO_low_depth_cutoff = mean_auto - sensitivity*sd_auto} 
      if(is.null(AUTO_high_depth_cutoff)) {AUTO_high_depth_cutoff  = mean_auto +sensitivity*sd_auto}
      XY_low_depth_cutoff = AUTO_low_depth_cutoff #Just define the XY cutoffs as the overall cutoffs
      XY_high_depth_cutoff = AUTO_high_depth_cutoff
    }
    cat("Mean depth cut-offs set for autosomes: mean depth must be >",AUTO_low_depth_cutoff,"and <",AUTO_high_depth_cutoff,"\n")
    cat("Mean depth cut-offs set for XY chromosomes: mean depth must be >",XY_low_depth_cutoff,"and <",XY_high_depth_cutoff,"\n")
  }
  filter_pass$mean_depth = sapply(1:nrow(COMB_mats$NV), assess_mean_depth, mat_list = COMB_mats, AUTO_low_depth_cutoff = AUTO_low_depth_cutoff, AUTO_high_depth_cutoff = AUTO_high_depth_cutoff, XY_low_depth_cutoff = XY_low_depth_cutoff, XY_high_depth_cutoff=XY_high_depth_cutoff)
  if(!is.null(pval_cutoff_dp2)) {filter_pass$pval_within_pos = ifelse(filter_params$pval_within_pos > pval_cutoff_dp2, 1, 0)} else {filter_pass$pval_within_pos = 1}
  if(!is.null(pval_cutoff_dp3)) {filter_pass$pval_within_pos_dp3 = ifelse(filter_params$pval_within_pos_dp3 > pval_cutoff_dp3, 1, 0)} else {filter_pass$pval_within_pos_dp3 = 1}
  filter_pass$max_depth_in_pos_samples = sapply(1:nrow(COMB_mats$NV), assess_max_depth_in_pos, mat_list = COMB_mats, min_depth_auto = min_depth_auto, min_depth_xy = min_depth_xy)
  filter_pass$max_pval_in_pos_sample = ifelse(filter_params$max_pval_in_pos_sample > min_pval_for_true_somatic, 1, 0)
  
  #Select the mutations for output
  select_muts = apply(filter_pass,1, function(x) all(x == 1))|rownames(filter_pass) %in% retain_muts
  filter_code = apply(filter_pass, 1, paste, collapse = "-")
  COMB_mats.tree.build = list_subset(COMB_mats, select_vector = select_muts)
  
  #Set the rownames to the mut_ref, and colnames to the sample names (without MTR or DEP)
  rownames(COMB_mats.tree.build$NV) = rownames(COMB_mats.tree.build$NR)  = rownames(COMB_mats.tree.build$PVal) <- COMB_mats.tree.build$mat$mut_ref
  colnames = gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats.tree.build$NV))
  colnames(COMB_mats.tree.build$NV) = colnames(COMB_mats.tree.build$NR) = colnames(COMB_mats.tree.build$PVal) <- colnames
  
  #BUILD THE GENOTYPE MATRIX
  #Select the "definite positive" samples by setting genotype to 1
  COMB_mats.tree.build$Genotype_bin = (COMB_mats.tree.build$NV >= min_variant_reads_SHARED & COMB_mats.tree.build$PVal > min_pval_for_true_somatic_SHARED) * 1
  #Select the "not sure" samples by setting genotype to 0.5
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV > 0 & COMB_mats.tree.build$PVal > 0.01 & COMB_mats.tree.build$Genotype_bin != 1] <- 0.5
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV >= 3 & COMB_mats.tree.build$PVal > 0.001 & COMB_mats.tree.build$Genotype_bin != 1] <- 0.5  
  COMB_mats.tree.build$Genotype_bin[(COMB_mats.tree.build$NV == 0) & (COMB_mats.tree.build$PVal > 0.1)] <- 0.5
  
  Genotype_shared_bin = COMB_mats.tree.build$Genotype_bin[rowSums(COMB_mats.tree.build$Genotype_bin == 1) > 1,]
  
  #Save the input parameters to a list
  params = list(input_set_ID = input_set_ID,
  				gender = gender,
  				retain_muts = retain_muts,
  				germline_pval_cutoff = -germline_pval_cutoff,
                rho_cutoff = rho_cutoff,
                XY_low_depth_cutoff = XY_low_depth_cutoff,
                XY_high_depth_cutoff = XY_high_depth_cutoff,
                AUTO_low_depth_cutoff= AUTO_low_depth_cutoff,
                AUTO_high_depth_cutoff= AUTO_high_depth_cutoff,
                pval_cutoff_dp2= pval_cutoff_dp2,
                pval_cutoff_dp3= pval_cutoff_dp3,
                min_depth_auto = min_depth_auto,
                min_depth_xy = min_depth_xy,
                min_pval_for_true_somatic = min_pval_for_true_somatic,
                min_variant_reads_SHARED = min_variant_reads_SHARED,
                min_pval_for_true_somatic_SHARED = min_pval_for_true_somatic_SHARED)
                
  #Save the summary stats of the run
  summary = data.frame(total_mutations = sum(select_muts),
  				total_SNVs = sum(COMB_mats.tree.build$mat$Mut_type == "SNV"),
  				total_INDELs = sum(COMB_mats.tree.build$mat$Mut_type == "INDEL"),
  				shared_mutations = nrow(Genotype_shared_bin),
  				shared_SNVs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "SNV"),
  				shared_INDELs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "INDEL"))
  #Extract the dummy dna_strings for tree building with MPBoot
  dna_strings = dna_strings_from_genotype(Genotype_shared_bin)
  
  #Print the run stats to the screen
  cat(summary$total_mutations,"total mutations\n")
  cat(summary$total_SNVs,"total SNVs\n")
  cat(summary$total_INDELs,"total INDELs\n")
  cat(summary$shared_mutations,"shared mutations\n")
  cat(summary$shared_SNVs,"shared SNVs\n")
  cat(summary$shared_INDELs,"shared INDELs\n")
  
  output = list(COMB_mats.tree.build = COMB_mats.tree.build, Genotype_shared_bin= Genotype_shared_bin, filter_code=filter_code, params = params, summary = summary, dna_strings = dna_strings)
  return(output)
}

#Function to create the dummy dna strings from the shared genotype matrix
dna_strings_from_genotype = function(genotype_mat) {
  Muts = rownames(genotype_mat)
  Ref = rep("W", length(Muts)) #W = Wild-type
  Alt = rep("V", length(Muts)) #V = Variant
  
  dna_strings = list()
  dna_strings[1] = paste(Ref, sep = "", collapse = "")
  
  for (k in 1:ncol(genotype_mat)){
    Mutations = Ref
    Mutations[genotype_mat[,k]==0.5] = '?'
    Mutations[genotype_mat[,k]==1] = Alt[genotype_mat[,k]==1]
    dna_string = paste(Mutations,sep="",collapse="")
    dna_strings[k+1]=dna_string
  }
  names(dna_strings)=c("Ancestral",colnames(genotype_mat))
  return(dna_strings)
}

