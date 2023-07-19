# select out a group of spots from infercnv dendrogram
library(infercnv)
library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)
library(hdf5r)
library(devtools)
library(SpatialInferCNV)

# provide inferncv directory, cluster name, and list of nodes belonging to this cluster
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
cluster_name <- args[2]

Consensus_AllBenigns <- read.dendrogram(file = paste0(dir, "infercnv.observations_dendrogram.txt"))

Consensus_AllBenigns_phylo <- as.phylo(Consensus_AllBenigns)

my.subtrees = subtrees(Consensus_AllBenigns_phylo)  

node_nums <- args[-c(1, 2)]
Merged <- data.frame(matrix(nrow = 0, ncol = 2))
for(i in 1:length(node_nums)){
  cur_node <- SelectingSubTreeData(my.subtrees, strtoi(node_nums[i]))
  Merged <- rbind(Merged, cur_node)
}

Merged$Node <- cluster_name
names(Merged)[2] <- "Histology"

write.csv(Merged, paste0(dir, "/", cluster_name, ".csv"), row.names = FALSE)
