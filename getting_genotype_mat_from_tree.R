details=filtered_muts$COMB_mats.tree.build$mat

#Copy the original genotype matrix to a "new_genotype_mat" that you will fill up based on the tree
new_genotype_mat = filtered_muts$COMB_mats.tree.build$Genotype_bin

#Set all the genotypes to 0 as a baseline
new_genotype_mat[,] <- 0

#Now go through each of the nodes in turn and replace the columns with 1's
for(i in tree$edge[,2]) {
  node=i
  info=get_edge_info(tree,details,node)
  samples=info$samples
  muts=details$mut_ref[info$idx]
  new_genotype_mat[muts,samples] <- 1
}


