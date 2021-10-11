nodeHeights = function (tree, ...) 
{
    if (hasArg(root.edge)) 
        root.edge <- list(...)$root.edge
    else root.edge <- FALSE
    if (root.edge) 
        ROOT <- if (!is.null(tree$root.edge)) 
            tree$root.edge
        else 0
    else ROOT <- 0
    nHeight <- function(tree) {
        tree <- reorder(tree)
        edge <- tree$edge
        el <- tree$edge.length
        res <- numeric(max(tree$edge))
        for (i in seq_len(nrow(edge))) res[edge[i, 2]] <- res[edge[i, 
            1]] + el[i]
        res
    }
    nh <- nHeight(tree)
    return(matrix(nh[tree$edge], ncol = 2L) + ROOT)
}

nodeheight=function (tree, node, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (node == (Ntip(tree) + 1)) 
    h <- 0
  else {
    a <- setdiff(c(getAncestors(tree, node), node), Ntip(tree) + 
                   1)
    h <- sum(tree$edge.length[sapply(a, function(x, e) which(e == 
                                                               x), e = tree$edge[, 2])])
  }
  h + ROOT
}

getAncestors=function (tree, node, type = c("all", "parent")) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  type <- type[1]
  if (type == "all") {
    aa <- vector()
    rt <- Ntip(tree) + 1
    currnode <- node
    while (currnode != rt) {
      currnode <- getAncestors(tree, currnode, "parent")
      aa <- c(aa, currnode)
    }
    return(aa)
  }
  else if (type == "parent") {
    aa <- tree$edge[which(tree$edge[, 2] == node), 1]
    return(aa)
  }
  else stop("do not recognize type")
}

getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

#Note - only apply the "get_edge_from_tree" option if starting from an SNV only tree. Function will assume that all the existing edge length is SNVs.
correct_edge_length = function(node, tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE, get_edge_from_tree=FALSE) {
  daughters <- getTips(tree = tree, node = node)
  #correct SNVs on edge, or set to 0 if want an indel only tree
  if(include_SNVs == TRUE) {
    if(get_edge_from_tree) {
    	nSNV=tree$edge.length[tree$edge[,2]==node]
    	} else {
    		nSNV = sum(details$node == node & details$Mut_type == "SNV")
    		}   		
    all_sens_SNVs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"SNV_sensitivity"]
    branch_SNV_sens = 1 - prod(1-all_sens_SNVs)  
    new_nSNV = nSNV/branch_SNV_sens
  } else {
    new_nSNV <- 0
    }
  #correct INDELs on edge, or set to 0 if want an SNV only tree
  if(include_indels == TRUE) {
    nINDEL = sum(details$node == node & details$Mut_type == "INDEL")
    all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
    branch_INDEL_sens = 1 - prod(1-all_sens_INDELs)
    new_nINDEL = nINDEL/branch_INDEL_sens
  } else {
    new_nINDEL <- 0
  }
  new_edge_length = new_nSNV + new_nINDEL
  return(new_edge_length)
}

get_subset_tree = function(tree, details, v.field = "Mut_type", value = "SNV") {
  get_new_edge_length = function(node, tree, details,v.field,value) {
    sum(details$node == node & details[v.field] == value)
  }
  tree_subset = tree
  tree_subset$edge.length = sapply(tree$edge[,2], get_new_edge_length, tree = tree, details = details,v.field = v.field,value=value)
  return(tree_subset)
}

get_corrected_tree = function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

get_mut_burden = function(tree) {
  mut_burden = nodeHeights(tree)[tree$edge[,2] %in% 1:length(tree$tip.label),2]
  return(mut_burden)
}

get_mut_burden_stats = function(tree) {
  mut_burden = get_mut_burden(tree)
  cat(paste("Mean mutation burden is", round(mean(mut_burden),digits = 1),"\n"))
  cat(paste("Range of mutation burden is", round(range(mut_burden)[1],digits = 1),"to",round(range(mut_burden)[2],digits = 1),"\n"))
  cat(paste("Standard deviation of mutation burden is", round(sd(mut_burden),digits = 1),"\n"))
}

#Function to calculate the absolute minimum number of clones by counting the number of times a parent node is shared, but
#a daughter node is recipient only (this would give you the number of extant transplanted clones if had full phylogeny)
get_minimum_clones=function(tree,donor_ID,recip_ID){
  shared_node_test=function(tree,node,donor_ID,recip_ID) {
    node_samples=getTips(tree,node)
    n_donor=sum(grepl(donor_ID,node_samples))
    n_recip=sum(grepl(recip_ID,node_samples))
    sharing_info=ifelse(n_donor>0&n_recip>0,"shared",ifelse(n_donor>0,"donor","recipient")) 
  }
  N=dim(tree$edge)[1]
  by_node=sapply(1:N,function(i) {
    node=tree$edge[i,2]
    sharing_info=shared_node_test(tree,node,donor_ID,recip_ID)
    if(sharing_info=="shared") {
      daughter_nodes=get_node_children(node,tree = tree)
      evidence_of_clone=0
      for(i in daughter_nodes) {
        daughter_sharing=shared_node_test(tree,node=i,donor_ID,recip_ID)
        if(daughter_sharing=="recipient") {
          evidence_of_clone=sum(1,evidence_of_clone)
        }
      }
    } else {
      evidence_of_clone=0
    }
    return(evidence_of_clone)
  })
  total_clones=sum(unlist(by_node))
  return(total_clones)
}

#Setup functions for AMOVA
# function to perform amova.
amova.fn <- function(distmat, groupnames, cell_key) {
  groupnums <- length(groupnames)
  dw <- c()
  cellnums <- c()
  for (i in 1:groupnums) {
    tgroup <- groupnames[i]
    tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type==tgroup])
    cellnums <- c(cellnums, length(tcells))
    dw <- c(dw, sum(distmat[tcells, tcells]))
  }
  tdist <- distmat[colnames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames], rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames]]
  dap <- (sum(tdist) - sum(dw))/2
  
  dfAP <- length(groupnames) - 1 
  dfWP <- sum(cellnums-1)
  N <- sum(cellnums)
  SSwp <- sum(dw/(cellnums*2))
  SSap <- sum(((dw + dap)/(2*N)) - (dw/(2*cellnums)))
  
  MSwp <- SSwp/dfWP
  MSap <- SSap/dfAP
  nc <- (N - (sum(cellnums^2)/N))/dfAP
  varwp <- MSwp
  varap <- (MSap - MSwp)/nc
  obsphi <- varap/(varwp + varap)
  return(obsphi)
}

# function to randomise sample labels and repeat
randamova.fn <- function(distmat, groupnames, cell_key) {
  
  # change added 2018.01.22: when randomizing, only include the part of the distance matrix that involves the cell types being considered.
  tcells <- which(rownames(distmat) %in% cell_key$Sample[cell_key$Cell_type %in% groupnames])
  #
  
  randmat <- distmat[tcells, tcells]
  colnames(randmat) <- sample(colnames(randmat))
  rownames(randmat) <- colnames(randmat)
  randphi <- amova.fn(distmat=randmat, groupnames=groupnames, cell_key=cell_key)
  return(randphi)
}

# function tying it all in together
amovapval.fn <- function(distmat, groupnames, cell_key, iterations, plottitle) {
  # calculate observed
  obsphi <- amova.fn(distmat=distmat, groupnames=groupnames, cell_key=cell_key)
  # calculate null
  randphis <- sapply(1:iterations, function(cell) randamova.fn(distmat = distmat, groupnames=groupnames, cell_key=cell_key))
  # calculate pval
  pval <- length(which(randphis>obsphi))/length(randphis) 
  
  hist(randphis, col="grey", 100, main=plottitle, xlab="Phi statistic")
  abline(v=obsphi, col="red", lwd=2)
  legend("topright", legend=paste0("Observed\n p = ", signif(pval,digits=2)), lwd=2, col="red", bty="n")
}


