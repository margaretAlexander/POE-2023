
format_spatial_data <- function(spatial_dir, sample_name) {
  counts <- ImportCountData(sample_name, paste0(sample_name, "/filtered_feature_bc_matrix.h5"))
  saveRDS(counts, file = paste0(sample_name, "/counts.RDS"))
  histology <- ImportHistologicalAnnotations(sample_name, paste0(sample_name, "/", histology_file)) 
  
  joined_counts <- MergingCountAndAnnotationData(sample_name, histology, counts)
  joined_counts <- joined_counts %>% column_to_rownames(var = "Genes")
  finalannot <- FinalAnnotations(histology, joined_counts)
  
  write.table(joined_counts, paste0(sample_name, "/joined_counts.tsv"), sep = "\t")
  write.table(finalannot, paste0(sample_name, "/final_annotations.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)
}

get_gene_ids <- function(){
  ensembl.ids <- rownames(joined_counts)
  gene.symbols <- gconvert(ensembl.ids, organism = "mmusculus", target = "MGI", filter_na = F)$target
  rem.ind <- which(duplicated(gene.symbols) == TRUE | is.na(gene.symbols))
  joined_counts <- joined_counts[-rem.ind, ]
  gene.symbols <- gene.symbols[-rem.ind]
  rownames(joined_counts) <- toupper(gene.symbols)
  write.table(joined_counts, paste0(sample_name, "/joined_counts2.tsv"), sep = "\t")
}

make_dendrogram <- function(infercnv_dir){
  all_obs <- read.dendrogram(file = paste0(infercnv_dir, "/infercnv.observations_dendrogram.txt"))
  
  all_obs_phylo <- as.phylo(all_obs)
  
  png(paste0(infercnv_dir, "/phylo_Nodes.png"),width=10000,height=2500, res = 300)
  plot(all_obs_phylo,show.tip.label = FALSE)
  nodelabels(text=1:all_obs_phylo$Nnode,node=1:all_obs_phylo$Nnode+Ntip(all_obs_phylo))
  dev.off()
}

run_sicnv <- function(spatial_dir, 
                      data_dir, 
                      reference = NULL, 
                      num_ref_groups = 0,
                      leiden_res = 0.05,
                      subculster = F,
                      supervised = F) {
  
  format_spatial_data()
  get_gene_ids()
                                                                                      
  if(supervised == F){
    
  }
  if(reference == F){
    
  }
  sicnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = paste0(data_dir, "/joined_counts2.tsv"),
                                              gene_order_file = paste0(data_dir, "mouse_gene_order.txt"),
                                              annotations_file = paste0(data_dir, "/final_annotations.tsv"), 
                                              delim = "\t", 
                                              ref_group_names = ref_groups,
                                              chr_exclude = c("chrM"))
  
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
  
  make_dendrogram()
  
}

map_back <- function(){
  all_obs <- read.dendrogram(file = paste0(dir, "/infercnv.observations_dendrogram.txt"))
  
  all_obs_phylo <- as.phylo(all_obs)
  
  my.subtrees = subtrees(all_obs_phylo)  
  
  node_nums <- args[-c(1, 2)]
  Merged <- data.frame(matrix(nrow = 0, ncol = 2))
  #colnames(Merged) <- c("Barcode", "Node")
  for(i in 1:length(node_nums)){
    cur_node <- SelectingSubTreeData(my.subtrees, strtoi(node_nums[i]))
    cur_node$Barcode <- sub("\\.", "-", cur_node$Barcode) 
    cur_node$Barcode <- sub(paste0(sample_name, "_"), "", cur_node$Barcode) # just don't add this in the original run...
    Merged <- rbind(Merged, cur_node)
  }
  
  names(Merged)[2] <- "Histology"
  
  write.csv(Merged, paste0(dir, "/infercnv_annotations.csv"), row.names = FALSE)
}