# plot spatially variable features and volcano plot 
library(Seurat)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(EnhancedVolcano)
library(dplyr)
library(fgsea)

# differential expression
args <- commandArgs(trailingOnly = TRUE)
spatial_dir <- args[1]
out_dir <- args[2]
# assuming spatial data is direct output from cellranger (ie "path/sample_name/outs")
sample <- strsplit(spatial_dir,"/")[[1]][length(strsplit(spatial_dir,"/")[[1]]) - 1]
seurat_obj <- Load10X_Spatial(spatial_dir)
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE) %>% 
  RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# can use sicnv clusters 
if(length(args) == 3){
  sicnv_ann <- read.csv(paste0("~/Documents/data/",sample, "/", args[3]))
  sicnv_ann$stroma[which(sicnv_ann$stroma == "")] <- "unlabeled"
  seurat_obj <- SetIdent(seurat_obj, sicnv_ann$Barcode, sicnv_ann$stroma)
}

# finding spatially variable features given annotations (either from seurat or sicnv)
fibro_hallmark_gene_sets <- read.csv("~/Documents/data/fibro_hallmark_gene_sets.applescript")
fibro_hallmark_gene_sets <- fibro_hallmark_gene_sets[which(fibro_hallmark_gene_sets$FDR < 0.05),]

for(cluster in levels(seurat_obj@active.ident)){
  markers <- FindMarkers(seurat_obj, test.use = "wilcox", ident.1 = cluster)
  write.csv(markers, paste0(out_dir, "/", sample, "/", cluster, "_markers.csv"))
  svg(paste0(out_dir, "/", sample, "/", cluster, "_spatial_feature_plot.svg"))
  print(SpatialFeaturePlot(seurat_obj, rownames(markers)[1:6], alpha = c(0.1, 1), ncol = 2))
  dev.off()
  svg(paste0(out_dir, "/", sample, "/", cluster, "_volcano.svg"), width = 12)
  print(EnhancedVolcano(markers,
                        lab = rownames(markers),
                        x = "avg_log2FC",
                        y = "p_val_adj",
                        title = paste0(sample, " ", cluster),
                        drawConnectors = TRUE,
                        col=c('gray', '#5FD35D', '#8D5FD3', '#D35F5F'),))
  dev.off()
}

