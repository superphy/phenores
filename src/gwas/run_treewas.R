# run_treewas.R: A script for running Treewas with different inputs
#
# Author: Matthew Whiteside
# Copyright 2018, Public Health Agency of Canada
# license: APL
# version: 2.0
# email: matthew.whiteside@canada.ca"

library(ape)
library(data.table)
library(treeWAS)

feature_file = snakemake@input[[1]]
pheno_file = snakemake@input[[2]]
tree_file = snakemake@input[[3]]
plot_file = snakemake@output[[2]]
r_file = snakemake@output[[1]]

loadTree <- function(treefile, do.root=FALSE, resolve.polytomies=TRUE) {
	# Load tree
	#
	# Root tree at midpoint.
	#
	# Args:
	#  treefile: Path to tree in newick format
	#  do.root: Boolean indicating if midpoint.root should be called on tree
	#  resolve.polytomies: Boolean indicating if multifurcating nodes should be expanded to bifurcating
	#
	# Returns:
	# 	phylo object
    
	otree = read.tree(treefile)
	tree = otree
	if(do.root) tree <- midpoint.root(tree)
	if(resolve.polytomies) {
		tree <- multi2di(tree,random=TRUE)
	}
	# Replace all zero length edges with a small value
	tree$edge.length[tree$edge.length < 1e-10] <- 1e-10

	tree
}

phenodf <- read.table(pheno_file, header=TRUE, sep='\t', as.is=TRUE)[c('sample','resistant')]
pheno <- phenodf[['resistant']]
names(pheno) = phenodf[['sample']]
featuredt <- fread(feature_file, header=TRUE, sep='\t')
featuremat <- as.matrix(featuredt[,-1,with=FALSE])
rownames(featuremat) <- featuredt[,1][[1]]
tree <- loadTree(tree_file)

out <- treeWAS(featuremat, pheno, tree, filename.plot = plot_file, mem.lim=snakemake@params[["mem"]])

saveRDS(out, file=r_file)




