# make_mash_groups.R: A script for computing clusters from a Mash distance matrix
#
# Author: Matthew Whiteside
# Copyright 2018, Public Health Agency of Canada
# license: APL
# version: 2.0
# email: matthew.whiteside@canada.ca"

library(FactoMineR)
library(dplyr)
library(readxl)

mash_dist_file = snakemake@input[[1]]
meta_excel_file = snakemake@input[[2]]

options(warn=-1)
metadf = read_excel(meta_excel_file, na='-')
options(warn=0)
distances <- read.csv(mash_dist_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)

# Identify population groups using hierarchical clustering on PCA-transformed mash distances
pcs <- PCA(distances, scale.unit=FALSE, ncp = 20, graph = FALSE)
clu <- HCPC(pcs, graph = FALSE, nb.clust=snakemake@params[["k"]])
clusterdf = clu$data.clust['clust']


# Break groups up into res/sus unless that would reduce group to < 10
samples = rownames(distances)
streptodf = data.frame(metadf[metadf$run %in% samples,c("run","phenotypic_streptomycin")], row.names="run")
colnames(streptodf) = c("str_res")
clusterdf = merge(clusterdf, streptodf, by=0, all=TRUE)
any(is.na(clusterdf))
clusterdf = clusterdf %>% add_count(clust, str_res)
clusterdf = clusterdf %>% mutate(stratgroup = case_when(
    n >= 10 ~ paste(clust,str_res,sep='-'), 
    n < 10 ~ ifelse(str_res == 0, paste(clust,1,sep='-'), paste(clust,0,sep='-'))
    ))

# Assign to one of 10 folds
set.seed(15)
clusterdf = clusterdf %>% sample_frac %>% group_by(stratgroup) %>% mutate(n = n(), nrow = row_number()-1)
clusterdf = clusterdf %>% mutate(SET = nrow %% 10) %>% ungroup

# Add filepath for future ease of use
clusterdf = clusterdf %>% mutate(fasta = paste('data/raw/genomes/',Row.names,'.fasta',sep=''))
clusterdf = clusterdf %>% as.data.frame

colnames(clusterdf) <- c('sample','cluster','resistant','group_n','group','group_row','fold','fasta')

write.csv(clusterdf, snakemake@output[[1]], quote=FALSE)


