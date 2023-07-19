# run infercnv on spatial trsncriptomic data and generate dendrogram
# use "stroma" from all samples to get larger reference since all samples are nearly entirely tumor
library(infercnv)
library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)
library(hdf5r)
library(devtools)
library(SpatialInferCNV)
library(gprofiler2)

args <- commandArgs(trailingOnly = TRUE)

# run from directory containing sample subdirectories 
obs_sample <- args[1] # sample that will make up the observation cells
options(digits = 10)
leiden_res_unsup <- as.double(args[2])
leiden_res_sup <- as.double(args[3])
sample_names <- list.dirs(getwd(), full.names = F, recursive = F) 
sample_names <- sample_names[sample_names != obs_sample]


# merging count data and histology
# histology file for observed sample should have the most "normal" spots labeled as "stroma", it can have other annotations
counts <- ImportCountData(obs_sample, paste0(obs_sample,"/filtered_feature_bc_matrix.h5"))
histology <- ImportHistologicalAnnotations(obs_sample, paste0(obs_sample, "/stroma_with_unlabeled.csv")) 
histology <- histology[which(histology$Histology != "stroma"),]
filtered_combo_ref <- read.csv("~/Documents/combo_ref_only_out/benign.csv")
histology <- merge(filtered_combo_ref, histology, all = T)
for(sample in sample_names){
  new_counts <- ImportCountData(sample, paste0(sample,"/filtered_feature_bc_matrix.h5")) 
  new_counts <- new_counts[which(new_counts$Barcode %in% filtered_combo_ref$Barcode),]
  if(nrow(new_counts) > 0){
    counts <- merge(new_counts, counts, all = T) 
  }
}
histology$Histology[which(histology$Histology == "")] <- "unknown"
saveRDS(counts, file = paste0(obs_sample, "/combo_counts.RDS"))

joined_counts <- MergingCountAndAnnotationData(paste0("merged_", obs_sample), histology, counts)
joined_counts <- joined_counts %>% column_to_rownames(var = "Genes")
finalannot <- FinalAnnotations(histology, joined_counts)

write.table(joined_counts, paste0(obs_sample, "/combo_joined_counts.tsv"), sep = "\t")
write.table(finalannot, paste0(obs_sample, "/combo_final_annotations.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)

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
write.table(joined_counts, paste0(obs_sample, "/combo_joined_counts2.tsv"), sep = "\t")


# #convert gene order file to ENSEMBL
# mouse_gene_order <- read.table(file = "data/mouse_gene_order.txt", stringsAsFactors = F)
# mouse_gene_order$V1 <- str_to_title(tolower(mouse_gene_order$V1))
# ensembl.rep <- gconvert(mouse_gene_order$V1, organism = "mmusculus", target = "ENSG", filter_na = F)$target


# make infercnv object
sicnv_combo_ref <- infercnv::CreateInfercnvObject(raw_counts_matrix = paste0(obs_sample, "/combo_joined_counts2.tsv"),
                                                  gene_order_file = "/Users/malexand/Documents/data/mouse_gene_order.txt",
                                                  annotations_file = paste0(obs_sample, "/combo_final_annotations.tsv"), 
                                                  delim = "\t", 
                                                  ref_group_names = "Normal",
                                                  chr_exclude = c("chrM"))


# unsupervised
sicnv_combo_ref_out_unsup <- infercnv::run(sicnv_combo_ref,
                                           cutoff = 0.1, # use 0.1 for 10x data
                                           out_dir = paste0(obs_sample, "/filtered_sicnv_combo_ref_unsup"),
                                           cluster_by_groups = F,
                                           HMM = F,
                                           denoise = T,
                                           analysis_mode = "subclusters",
                                           write_phylo = T,
                                           num_threads = 12,
                                           leiden_resolution = 0.001,
                                           resume_mode = F)

all_obs <- read.dendrogram(file = paste0(obs_sample, "/sicnv_combo_ref_unsup/infercnv.observations_dendrogram.txt"))

all_obs_phylo <- as.phylo(all_obs)

png(paste0(obs_sample, "/sicnv_combo_ref_unsup/phylo_Nodes.png"),width=10000,height=2500, res = 300)
plot(all_obs_phylo,show.tip.label = FALSE)
nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
dev.off()

# supervised
options(scipen = 1000)
sicnv_combo_ref_out_sup <- suppressWarnings(infercnv::run(sicnv_combo_ref,
                                                          cutoff = 0.1, # use 0.1 for 10x data
                                                          out_dir = paste0(obs_sample, "/sicnv_combo_ref_sup"),
                                                          cluster_by_groups = F,
                                                          HMM = T,
                                                          denoise = T,
                                                          analysis_mode = "subcluster", 
                                                          write_phylo = T,
                                                          num_threads = 12,
                                                          leiden_resolution = leiden_res_sup,
                                                          resume_mode = F))

all_obs <- read.dendrogram(file = paste0(obs_sample, "/sicnv_combo_ref_sup/infercnv.observations_dendrogram.txt"))

all_obs_phylo <- as.phylo(all_obs)

png(paste0(obs_sample, "/sicnv_combo_ref_sup/phylo_Nodes.png"),width=10000,height=2500, res = 300)
plot(all_obs_phylo,show.tip.label = FALSE)
nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
dev.off()

# examine reference
combo_ref <- histology[which(histology$Histology == "stroma"),]
ref_counts <- joined_counts[,which(colnames(joined_counts) %in% histology$Barcode)]
finalannot_ref_only <- FinalAnnotations(combo_ref, ref_counts)
write.table(finalannot_ref_only, "~/Documents/combo_ref_only.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ref_counts, "~/Documents/combo_ref_counts.tsv", sep = "\t")

sicnv_combo_ref_only <- infercnv::CreateInfercnvObject(raw_counts_matrix = "~/Documents/combo_ref_counts.tsv",
                                                       gene_order_file = "/Users/malexand/Documents/data/mouse_gene_order.txt",
                                                       annotations_file = "~/Documents/combo_ref_only.tsv", 
                                                       delim = "\t", 
                                                       ref_group_names = NULL,
                                                       chr_exclude = c("chrM"))

sicnv_combo_ref_only_out <- infercnv::run(sicnv_combo_ref_only,
                                          cutoff = 0.1, # use 0.1 for 10x data
                                          out_dir = "~/Documents/combo_ref_only_out",
                                          cluster_by_groups = F,
                                          HMM = F,
                                          denoise = T,
                                          analysis_mode = "samples",
                                          write_phylo = T,
                                          num_threads = 12,
                                          resume_mode = F)

all_obs <- read.dendrogram(file = "~/Documents/combo_ref_only_out/infercnv.observations_dendrogram.txt")

all_obs_phylo <- as.phylo(all_obs)

png("~/Documents/combo_ref_only_out/phylo_Nodes.png",width=10000,height=2500, res = 300)
plot(all_obs_phylo,show.tip.label = FALSE)
nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
dev.off()
