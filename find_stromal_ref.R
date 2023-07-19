# find most "normal" cells from stromal groups
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

get_gene_symbols <-function(joined_counts, data_dir){
  ensembl.ids <- rownames(joined_counts)
  gene.symbols <- gconvert(ensembl.ids, organism = "mmusculus", target = "MGI", filter_na = F)$target
  #rownames(joined_counts) <- gene.symbols #non unique
  #remove NAs
  rem.ind <- which(duplicated(gene.symbols) == TRUE | is.na(gene.symbols))
  joined_counts <- joined_counts[-rem.ind, ]
  gene.symbols <- gene.symbols[-rem.ind]
  rownames(joined_counts) <- toupper(gene.symbols)
  write.table(joined_counts, paste0(data_dir, "joined.counts2.tsv"), sep = "\t")
}

# accept path to histology file and output directory as command line argument

data_dir <- "./data/mixed_regions/"
B.counts <- ImportCountData("B", paste0(data.dir, "filtered_feature_bc_matrix.h5"))
saveRDS(B.counts, file = paste0(data_dir, "/B.counts.RDS"))

histology <- ImportHistologicalAnnotations("B", paste0(data_dir, "Histology_no_unlabeled.csv")) 
histology <- histology[which(histology$Histology %in% c("stromal_mixed1", "stromal_mixed2")),]
#histology$Histology[which(histology$Histology == "")] <- "unknown"

joined_counts <- MergingCountAndAnnotationData("B", histology, B.counts)
joined_counts <- joined_counts %>% column_to_rownames(var = "Genes")
finalannot <- FinalAnnotations(histology, joined_counts)

write.table(joined_counts, paste0(data_dir, "joined.counts.tsv"), sep = "\t")
write.table(finalannot, paste0(data_dir, "final_annotations.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)

# library(biomaRt)
# mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
# genes <- colnames(B.counts)
# genes <- genes[-1]
# glist <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "entrezgene_id", "description", "external_gene_name"), values = genes, mart = mart)
# glist$external_gene_name <- toupper(glist$external_gene_name)
# rownames(glist) <- glist$ensembl_gene_id

get_gene_symbols(joined_counts, data_dir)


# #convert gene order file to ENSEMBL
# mouse_gene_order <- read.table(file = "data/mouse_gene_order.txt", stringsAsFactors = F)
# mouse_gene_order$V1 <- str_to_title(tolower(mouse_gene_order$V1))
# ensembl.rep <- gconvert(mouse_gene_order$V1, organism = "mmusculus", target = "ENSG", filter_na = F)$target



infercnv_mixedregions_stromal <- infercnv::CreateInfercnvObject(raw_counts_matrix = paste0(data_dir, "joined.counts2.tsv"),
                                                        gene_order_file = "data/mouse_gene_order.txt",
                                                        annotations_file = paste0(data_dir, "final_annotations.tsv"), 
                                                        delim = "\t", 
                                                        ref_group_names = NULL,
                                                        chr_exclude = c("chrM"))


infercnv_mixedregions_stromal_out <- infercnv::run(infercnv_mixedregions_stromal,
                                           cutoff = 0.1,
                                           out_dir = "infercnv_out/mixed_regions/with_unlabeled/find_stromal_ref",
                                           cluster_by_groups = F,
                                           HMM = F,
                                           denoise = T,
                                           analysis_mode = "samples",
                                           write_phylo = T,
                                           num_threads = 12)

Consensus_AllBenigns <- read.dendrogram(file = "./infercnv_out/mixed_regions/with_unlabeled/find_stromal_ref/infercnv.observations_dendrogram.txt")

Consensus_AllBenigns_phylo <- as.phylo(Consensus_AllBenigns)

my.subtrees = subtrees(Consensus_AllBenigns_phylo)  

png("stromal_phylo_Nodes.png",width=10000,height=2500, res = 300)
plot(Consensus_AllBenigns_phylo,show.tip.label = FALSE)
nodelabels(text=1:Consensus_AllBenigns_phylo$Nnode,node=1:Consensus_AllBenigns_phylo$Nnode+Ntip(Consensus_AllBenigns_phylo))
dev.off()

purest_benign <- read.csv( "./stromal_phylo_Nodes.csv")

joined_counts <- MergingCountAndAnnotationData("B", purest_benign, B.counts)
joined_counts <- joined_counts %>% column_to_rownames(var = "Genes")
finalannot <- FinalAnnotations(histology, joined_counts)

write.table(joined_counts, paste0(data_dir, "joined.counts.tsv"), sep = "\t")
write.table(finalannot, paste0(data_dir, "final_annotations.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)

get_gene_symbols(joined_counts, data_dir)

# check that no significant cnv signals are detected with new refercne group

infercnv_mixedregions <- infercnv::CreateInfercnvObject(raw_counts_matrix = paste0(data_dir, "joined.counts2.tsv"),
                                                        gene_order_file = "data/mouse_gene_order.txt",
                                                        annotations_file = paste0(data_dir, "final_annotations.tsv"), 
                                                        delim = "\t", 
                                                        ref_group_names = NULL,
                                                        chr_exclude = c("chrM"))

infercnv_mixedregions_out <- infercnv::run(infercnv_mixedregions,
                                           cutoff = 0.1,
                                           out_dir = "infercnv_out/mixed_regions/with_unlabeled/test_stromal_ref",
                                           cluster_by_groups = F,
                                           HMM = F,
                                           denoise = T,
                                           analysis_mode = "samples",
                                           write_phylo = T,
                                           num_threads = 12)
