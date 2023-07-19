library(Seurat)
library(biomaRt)
library(gprofiler2)
library(biomaRt)
library(readxl)

tpm <- function(raw_counts, gene_lengths) {
  
  rpk <- raw_counts*1e3 / gene_lengths
  pm <- colSums(rpk) / 1e6
  return(log(t(t(rpk) / pm) + 1, 2))
  
}

# makes scatter plot with spatial vs single cell data for given sample
args <- commandArgs(trailingOnly = TRUE)
spatial_data <- args[1]
out_dir <- args[2]

mes_and_adrn_markers <- read_excel("~/Documents/data/mes_and_adrn_markers.xlsx")[-1,]
colnames(mes_and_adrn_markers) <- c("gene", "marker")
markers <- data.frame(mes_and_adrn_markers$marker, row.names = mes_and_adrn_markers$gene)

sample <- strsplit(spatial_data,"/")[[1]][length(strsplit(spatial_data,"/")[[1]]) - 1]
seurat_obj <- Load10X_Spatial(spatial_data)
st_counts <- seurat_obj@assays[["Spatial"]]@counts
sicnv_ann <- read.csv("~/Documents/data/NB849_B/stroma_with_unlabeled.csv")
sicnv_ann$stroma[which(sicnv_ann$stroma == "")] <- "unlabeled"

gene_symbols <- toupper(st_counts@Dimnames[[1]])
ensembl_ids <- gconvert(gene_symbols, organism = "mmusculus", target = "ENSG", filter_na = F)$target
mouse <- useEnsembl("genes", "mmusculus_gene_ensembl")
gene_coords <- getBM(attributes=c("mgi_symbol","ensembl_gene_id", "start_position","end_position"), filters="ensembl_gene_id", values=ensembl_ids, mart=mouse)
gene_coords$size <- gene_coords$end_position - gene_coords$start_position
gene_coords$mgi_symbol <- toupper(gene_coords$mgi_symbol)
gene_coords <- gene_coords[order(gene_coords$mgi_symbol),]
gene_coords <- gene_coords[which(gene_coords$mgi_symbol %in% gene_symbols),]
gene_coords <- gene_coords[which(!(duplicated(gene_coords$mgi_symbol))),]

st_counts <- as.data.frame(as.matrix(st_counts))
st_tpm <- tpm(st_counts, gene_coords$size)

markers_boxplot <- function(st_tpm, sicnv_ann, cluster, markers){
  cluster_tpm <- st_tpm[, which(colnames(st_tpm) %in% sicnv_ann$Barcode[which(sicnv_ann$stroma == cluster)])]
  cluster_tpm_avg <- data.frame(rowMeans(cluster_tpm))
  rownames(cluster_tpm_avg) <- toupper(rownames(cluster_tpm_avg))
  
  marker_tpm <- merge(cluster_tpm_avg, markers, by = 0, all = T)
  marker_tpm$mes_and_adrn_markers.marker[which(is.na(marker_tpm$mes_and_adrn_markers.marker))] <- "Neither"
  marker_tpm <- marker_tpm[-which(is.na(marker_tpm$rowMeans.cluster_tpm.)),]
  marker_tpm$mes_and_adrn_markers.marker <- as.factor(marker_tpm$mes_and_adrn_markers.marker)
  
  boxplot(rowMeans.cluster_tpm.~mes_and_adrn_markers.marker, marker_tpm, xlab = NULL, ylab = NULL)
}
stroma_tpm <- st_tpm[, which(colnames(st_tpm) %in% sicnv_ann$Barcode[which(sicnv_ann$stroma == "stroma")])]
unlabeled_tpm <- st_tpm[, which(colnames(st_tpm) %in% sicnv_ann$Barcode[which(sicnv_ann$stroma == "unlabeled")])]

stroma_tpm_avg <- data.frame(rowMeans(stroma_tpm))
rownames(stroma_tpm_avg) <- toupper(rownames(stroma_tpm_avg))



unlabeled_tpm_avg <- data.frame(rowMeans(unlabeled_tpm))
rownames(unlabeled_tpm_avg) <- toupper(rownames(unlabeled_tpm_avg))

marker_tpm <- merge(st_tpm_avg, markers, by = 0, all = T)
marker_tpm$mes_and_adrn_markers.marker[which(is.na(marker_tpm$mes_and_adrn_markers.marker))] <- "Neither"
marker_tpm <- marker_tpm[-which(is.na(marker_tpm$rowMeans.st_tpm.)),]
marker_tpm$mes_and_adrn_markers.marker <- as.factor(marker_tpm$mes_and_adrn_markers.marker)

boxplot(rowMeans.st_tpm.~mes_and_adrn_markers.marker, marker_tpm, xlab = NULL, ylab = NULL)

for(cluster in unique(sicnv_ann$stroma)){
  markers_boxplot(st_tpm, sicnv_ann, cluster, markers)
}

# boxplots comparing stroma and unlabeled
adrn <- mes_and_adrn_markers$gene[which(mes_and_adrn_markers$marker == "ADRN")]
mes <- mes_and_adrn_markers$gene[which(mes_and_adrn_markers$marker == "MES")]

stroma_mes <- data.frame(stroma_tpm_avg[which(rownames(stroma_tpm_avg) %in% mes),])
colnames(stroma_mes) <- "avg_log_tpm"
stroma_mes$label <- "Normal"
stroma_adrn <- stroma_tpm_avg[which(rownames(stroma_tpm_avg) %in% adrn),]
stroma_neither <- stroma_tpm_avg[which(!(rownames(stroma_tpm_avg) %in% mes | rownames(stroma_tpm_avg) %in% adrn)),]

unlabeled_mes <- data.frame(unlabeled_tpm_avg[which(rownames(unlabeled_tpm_avg) %in% mes),])
colnames(unlabeled_mes) <- "avg_log_tpm"
unlabeled_mes$label <- "Non-Normal"
unlabeled_adrn <- unlabeled_tpm_avg[which(rownames(unlabeled_tpm_avg) %in% adrn),]
unlabeled_neither <- unlabeled_tpm_avg[which(!(rownames(unlabeled_tpm_avg) %in% mes | rownames(unlabeled_tpm_avg) %in% adrn)),]

mes_tpm <- merge(stroma_mes, unlabeled_mes, all = T)
mes_tpm$label <- as.factor(mes_tpm$label)
boxplot(avg_log_tpm~label, mes_tpm, xlab = NULL, ylab = NULL, widths = c(0.1, 0.1) )
# fisher test

stroma_markers <- read.csv("~/Documents/output/NB849_B/stroma_markers.csv")
stroma_markers <- stroma_markers[which(stroma_markers$p_val_adj < 0.05),]
stroma_up <- stroma_markers[which(stroma_markers$avg_log2FC >= 1 & stroma_markers$p_val_adj < 0.05),]
stroma_up <- toupper(stroma_up$X)
stroma_not_up <- toupper(stroma_markers$X[which(!(toupper(stroma_markers$X) %in% stroma_up))])

stroma_up_and_mes <- length(which(mes %in% stroma_up))
stroma_up_and_not_mes <- length(stroma_up) - stroma_up_and_mes

stroma_not_up_and_mes <- length(which(mes %in% stroma_not_up))
stroma_not_up_and_not_mes <- length(stroma_not_up) - stroma_not_up_and_mes

fisher.test(matrix(c(stroma_up_and_mes, stroma_up_and_not_mes, stroma_not_up_and_mes, stroma_not_up_and_not_mes), nrow=2))



