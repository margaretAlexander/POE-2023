# run infercnv on sc RNA-seq data for specific smaples in merged seurat object
library(infercnv)
library(dplyr)
library(Seurat)
library(SingleR)
library(celldex)


# command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
out_dir <- args[2]

# load seurat object
mouse_sc_rna_seq <- readRDS(paste0(input_dir, "/NB.DBHiCre.merged.RDS"))
counts <- mouse_sc_rna_seq@assays[["RNA"]]@counts
ref <- celldex::MouseRNAseqData()


sample_names <-  levels(mouse_sc_rna_seq@active.ident)
sample_nums <- gsub("NB", "", sample_names)

for(sample in sample_nums){
  sample_counts <- counts[,which(grepl(sample, colnames(counts)))]
  
  pred_ann <- SingleR(test = sample_counts, ref = ref, labels = ref$label.fine, assay.type.test=1)
  sample_ann <- data.frame(pred_ann$labels, row.names = rownames(pred_ann))
  
  # create infercnv object
  infercnvobj = CreateInfercnvObject(raw_counts_matrix = sample_counts,
                                     annotations_file = sample_ann,
                                     delim = "\t",
                                     gene_order_file = paste0(input_dir, "/mouse_gene_order.txt"),
                                     ref_group_names = c("Endothelial cells", "T cells", "B cells"))
  
  # run infercnv on previous object, may need to suppress warnings if running on subclusters mode
  infercnvres = infercnv::run(infercnvobj,
                              cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir = paste0(out_dir, "/infercnv_", sample), 
                              cluster_by_groups = F,  
                              denoise = T,
                              HMM = F,
                              num_threads = 12,
                              analysis_mode = "samples",
                              write_phylo = T) 
}

