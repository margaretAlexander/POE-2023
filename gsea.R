# Run gsea on spatial data
library(escape)
library(dittoSeq)
library(msigdbr)
library(Seurat)

tpm <- function(raw_counts, gene_lengths) {
  
  rpk <- raw_counts*1e3 / gene_lengths
  pm <- colSums(rpk) / 1e6
  return(t(t(rpk) / pm) + 1)
  
}


hallmark_pathwaysDF <- msigdbr("Mus musculus", category = "H")
hallmark_pathways <- split(as.character(hallmark_pathwaysDF$gene_symbol), hallmark_pathwaysDF$gs_name)

kegg_pathwaysDF <- msigdbr("Mus musculus", category = "C2", subcategory = "CP:KEGG")
kegg_pathways <- split(as.character(kegg_pathwaysDF$gene_symbol), kegg_pathwaysDF$gs_name)

args <- commandArgs(trailingOnly = TRUE)
spatial_dir<- args[1]
out_dir <- args[2]
# Assuming spatial data is direct output from cellranger (ie "path/sample_name/outs")
sample <- strsplit(spatial_dir,"/")[[1]][length(strsplit(spatial_dir,"/")[[1]]) - 1]
seurat_obj <- Load10X_Spatial(spatial_dir)

if(length(args) == 3){
  sicnv_ann <- read.csv(paste0("~/Documents/data/",sample, "/", args[3]))
  sicnv_ann$stroma[which(sicnv_ann$stroma == "")] <- "unlabeled"
  colnames(sicnv_ann) <- c("barcodes", "Histology")
  seurat_obj <- SetIdent(seurat_obj, sicnv_ann$barcodes, sicnv_ann$Histology)
  seurat_obj<- Seurat::AddMetaData(seurat_obj, Idents(seurat_obj), col.name = "annotations")
}
seurat_obj <- RenameAssays(seurat_obj, Spatial = "RNA") 

hallmark_seurat <- enrichIt(obj = seurat_obj, 
                      gene.sets = hallmark_pathways, 
                      groups = 1000, cores = 12, 
                      min.size = 5)
seurat_obj<- Seurat::AddMetaData(seurat_obj, hallmark_seurat)

kegg_seurat <- enrichIt(obj = seurat_obj, 
                        gene.sets = kegg_pathways, 
                        groups = 1000, cores = 12, 
                        min.size = 5)
seurat_obj<- Seurat::AddMetaData(seurat_obj, kegg_seurat)

colors <- colorRampPalette(c("#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF"))

# Fibroblast hallmark gene sets obtained from pscb.stjude.org
fibro_hallmark_gene_sets <- read.csv("~/Documents/data/fibro_hallmark_gene_sets.applescript")
fibro_hallmark_gene_sets <- fibro_hallmark_gene_sets[which(fibro_hallmark_gene_sets$FDR < 0.05),]

#Heat map
dittoHeatmap(seurat_obj, genes = NULL, metas = fibro_hallmark_gene_sets$X, 
             annot.by = "annotations", 
             fontsize = 7, 
             heatmap.colors = colors(50))

# Boxplot
multi_dittoPlot(seurat_obj, vars = fibro_hallmark_gene_sets$X, 
                group.by = "annotations", plots = c("boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

