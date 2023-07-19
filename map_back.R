# Based on method from Erickson et al.

# Map sicnv identified groups back to slide
# Have to manually determine which nodes to use from dendrogram and use nodes numbers as parameters

library(biomaRt)
library(infercnv)
library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)
library(hdf5r)
library(devtools)
library(SpatialInferCNV)
library(gprofiler2)

# args should be infercnv directory followed by nodes of subclusters
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
sample_name <- args[2]

all_obs <- read.dendrogram(file = paste0(dir, "/infercnv.observations_dendrogram.txt"))

all_obs_phylo <- as.phylo(all_obs)

my.subtrees = subtrees(all_obs_phylo)  

node_nums <- args[-c(1, 2)]
Merged <- data.frame(matrix(nrow = 0, ncol = 2))
for(i in 1:length(node_nums)){
  cur_node <- SelectingSubTreeData(my.subtrees, strtoi(node_nums[i]))
  cur_node$Barcode <- sub("\\.", "-", cur_node$Barcode) 
  cur_node$Barcode <- sub(paste0(sample_name, "_"), "", cur_node$Barcode) 
  Merged <- rbind(Merged, cur_node)
}

names(Merged)[2] <- "Histology"

write.csv(Merged, paste0(dir, "/infercnv_annotations.csv"), row.names = FALSE)
# Import annotations into loupe browser to look for morphological differences and spatial context of clusters
