# run infercnv on spatial trsncriptomic data and generates dendrogram, runs both with and without hmm
library(infercnv)
library(tidyverse)
library(Seurat)
library(phylogram)
library(ape)
library(hdf5r)
library(devtools)
library(SpatialInferCNV)
library(gprofiler2)

# command line arguments
# run from data directory
args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1] # sample name
output_dir <- args[2] # infercnv output directory
histology_file <- args[3]# file with cell type annotations/cell groupings (probably from loupe browser)
leiden_res_sup <- args[4]
leiden_res_unsup <- args[5]
ref_groups <- args[6:length(args)]

counts <- ImportCountData(sample_name, paste0(sample_name, "/filtered_feature_bc_matrix.h5"))
saveRDS(counts, file = paste0(sample_name, "/counts.RDS"))
histology <- ImportHistologicalAnnotations(sample_name, paste0(sample_name, "/", histology_file)) 

joined_counts <- MergingCountAndAnnotationData(sample_name, histology, counts)
joined_counts <- joined_counts %>% column_to_rownames(var = "Genes")
finalannot <- FinalAnnotations(histology, joined_counts)

write.table(joined_counts, paste0(sample_name, "/joined_counts.tsv"), sep = "\t")
write.table(finalannot, paste0(sample_name, "/final_annotations.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)

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
write.table(joined_counts, paste0(sample_name, "/joined_counts2.tsv"), sep = "\t")


# #convert gene order file to ENSEMBL
# mouse_gene_order <- read.table(file = "data/mouse_gene_order.txt", stringsAsFactors = F)
# mouse_gene_order$V1 <- str_to_title(tolower(mouse_gene_order$V1))
# ensembl.rep <- gconvert(mouse_gene_order$V1, organism = "mmusculus", target = "ENSG", filter_na = F)$target



sicnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = paste0(sample_name, "/joined_counts2.tsv"),
                                            gene_order_file = "mouse_gene_order.txt",
                                            annotations_file = paste0(sample_name, "/final_annotations.tsv"), 
                                            delim = "\t", 
                                            ref_group_names = ref_groups,
                                            chr_exclude = c("chrM")
                                                                )

sicnv_unsup <- infercnv::run(sicnv_obj,
                            cutoff = 0.1, # use 0.1 for 10x data
                            out_dir = paste0(output_dir, "/", sample_name, "/sicnv_unsup"),
                            cluster_by_groups = F,
                            HMM = F,
                            denoise = T,
                            analysis_mode = "subcluster", 
                            write_phylo = T,
                            num_threads = 12,
                            leiden_resolution = leiden_res_unsup,
                            num_ref_groups = length(ref_groups))

all_obs <- read.dendrogram(file = paste0(output_dir, "/", sample_name, "/sicnv_unsup/infercnv.observations_dendrogram.txt"))

all_obs_phylo <- as.phylo(all_obs)

png(paste0(output_dir, "/", sample_name, "/sicnv_unsup/phylo_Nodes.png"),width=10000,height=2500, res = 300)
plot(all_obs_phylo,show.tip.label = FALSE)
nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
dev.off()

sicnv_sup <- infercnv::run(sicnv_obj,
                           cutoff = 0.1, # use 0.1 for 10x data
                           out_dir = paste0(output_dir, "/", sample_name, "/sicnv_sup"),
                           cluster_by_groups = T,
                           HMM = T,
                           denoise = T,
                           analysis_mode = "subcluster", 
                           write_phylo = T,
                           num_threads = 12,
                           leiden_resolution = leiden_res_sup,
                           num_ref_groups = length(ref_groups))


all_obs <- read.dendrogram(file = paste0(output_dir, "/", sample_name, "/sicnv_sup/infercnv.observations_dendrogram.txt"))

all_obs_phylo <- as.phylo(all_obs)

png(paste0(output_dir, "/", sample_name, "/sicnv_sup/phylo_Nodes.png"),width=10000,height=2500, res = 300)
plot(all_obs_phylo,show.tip.label = FALSE)
nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
dev.off()

