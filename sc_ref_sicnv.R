# Using sc reference with st obs
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library(infercnv)
library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)
library(hdf5r)
library(devtools)
install_github("aerickso/SpatialInferCNV")
library(SpatialInferCNV)
library(gprofiler2)
library(dplyr)
library(SingleR)
library(celldex)

# get sc RNA-seq data
# run from data directory
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
sample_num <- args[2]
samples <- list.dirs(data_dir)
sample_name <- samples[which(grepl(sample_num, samples))]
leiden_res_sup <- args[3]
leiden_res_unsup <- args[4]

mouse_sc_rna_seq <- readRDS(paste0(data_dir, "NB.DBHiCre.merged.RDS"))
counts <- mouse_sc_rna_seq@assays[["RNA"]]@counts
ref <- celldex::MouseRNAseqData()

# separate out NB849
sc_counts <- as.data.frame(as.matrix(counts[,which(grepl(sample_num, colnames(counts)))]))

pred_ann <- SingleR(test = sc_counts, ref = ref, labels = ref$label.fine, assay.type.test=1)
sc_ann <- data.frame(rownames(pred_ann), pred_ann$labels)
colnames(sc_ann) <- c("Barcode", "Histology")
sc_ref <- sc_ann %>% filter(sc_ann$Histology %in% c("Endothelial cells", "T cells", "B cells"))
sc_counts <- sc_counts[, which(colnames(sc_counts) %in% sc_ref$Barcode)]

# get spatial transcriptomic data
histology_file <- "all_unlabeled.csv"
st_counts <- ImportCountData(sample_name, paste0(sample_name, "/filtered_feature_bc_matrix.h5"))
saveRDS(st_counts, file = paste0(sample_name, "/counts.RDS"))

histology <- ImportHistologicalAnnotations(sample_name, paste0(sample_name, "/", histology_file)) 
#histology$Histology[which(histology$Histology == "")] <- "unknown"

joined_counts <- MergingCountAndAnnotationData(sample_name, histology, st_counts)
joined_counts <- joined_counts %>% column_to_rownames(var = "Genes")

finalannot <- FinalAnnotations(histology, joined_counts)
finalannot <- merge(finalannot, sc_ref, all = T)

write.table(joined_counts, paste0(sample_name, "/joined.counts.sc.st.tsv"), sep = "\t")
write.table(finalannot, paste0(sample_name, "/final_annotations_sc_st.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)

# library(biomaRt)
# mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
# genes <- colnames(B.counts)
# genes <- genes[-1]
# glist <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "entrezgene_id", "description", "external_gene_name"), values = genes, mart = mart)
# glist$external_gene_name <- toupper(glist$external_gene_name)
# rownames(glist) <- glist$ensembl_gene_id

ensembl.ids <- rownames(joined_counts)
gene.symbols <- gconvert(ensembl.ids, organism = "mmusculus", target = "MGI", filter_na = F)$target
#rownames(joined_counts) <- gene.symbols #non unique
#remove NAs
rem.ind <- which(duplicated(gene.symbols) == TRUE | is.na(gene.symbols))
joined_counts <- joined_counts[-rem.ind, ]
gene.symbols <- gene.symbols[-rem.ind]
rownames(joined_counts) <- toupper(gene.symbols)

# filter out genes that don't occur in both assays
joined_counts <- joined_counts %>% filter(rownames(joined_counts) %in% rownames(sc_counts))
sc_counts <- sc_counts %>% filter(rownames(sc_counts) %in% rownames(joined_counts))

sc_st_counts <- merge(x = joined_counts, y = sc_counts, by = 0, all = T)
rownames(sc_st_counts) <- sc_st_counts$Row.names
sc_st_counts <- sc_st_counts[,-1]

write.table(sc_st_counts, paste0( sample_name, "/sc_st_counts.tsv"), sep = "\t")

sicnv_with_sc_ref <- infercnv::CreateInfercnvObject(raw_counts_matrix = paste0(sample_name, "/sc_st_counts.tsv"),
                                                                gene_order_file = "mouse_gene_order.txt",
                                                                annotations_file = paste0(sample_num, "/final_annotations_sc_st.tsv"), 
                                                                delim = "\t", 
                                                                ref_group_names = c("Endothelial cells", "T cells", "B cells"),
                                                                chr_exclude = c("chrM"))


sicnv_with_sc_ref_unsup<- infercnv::run(sicnv_with_sc_ref,
                                    cutoff = 0.1,
                                    out_dir = paste0(out_dir, "/", sample_name,"/sicnv_with_sc_ref_unsup"),
                                    cluster_by_groups = F,
                                    HMM = F,
                                    denoise = T,
                                    analysis_mode = "subcluster",
                                    write_phylo = T,
                                    num_threads = 12,
                                    leiden_resolution = leiden_res_unsup)

all_obs <- read.dendrogram(file = paste0(out_dir, "/", sample_name, "/sicnv_with_sc_ref_unsup/infercnv.observations_dendrogram.txt"))

all_obs_phylo <- as.phylo(all_obs)

png(paste0(out_dir, "/", sample_name, "/sicnv_with_sc_ref_unsup/phylo_Nodes.png"),width=10000,height=2500, res = 300)
plot(all_obs_phylo,show.tip.label = FALSE)
nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
dev.off()

sicnv_with_sc_ref_sup<- infercnv::run(sicnv_with_sc_ref,
                                        cutoff = 0.1,
                                        out_dir = paste0(out_dir, "/", sample_name,"/sicnv_with_sc_ref_sup"),
                                        cluster_by_groups = T,
                                        HMM = T,
                                        denoise = T,
                                        analysis_mode = "subcluster",
                                        write_phylo = T,
                                        num_threads = 12,
                                        leiden_resolution = leiden_res_sup)

all_obs <- read.dendrogram(file = paste0(out_dir, "/", sample_name, "/sicnv_with_sc_ref_sup/infercnv.observations_dendrogram.txt"))

all_obs_phylo <- as.phylo(all_obs)

png(paste0(out_dir, "/", sample_name, "/sicnv_with_sc_ref_sup/phylo_Nodes.png"),width=10000,height=2500, res = 300)
plot(all_obs_phylo,show.tip.label = FALSE)
nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
dev.off()



