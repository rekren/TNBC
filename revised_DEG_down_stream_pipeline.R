setwd("D:/work/tnbc_cell/mapped/svm")
countData <- read.table(file = "whole_exp.txt" , header = TRUE, sep = '\t',row.names = 1)
countData <- countData[,sort(names(countData))]
names(countData) <-(gsub("\\.", "\\-", names(countData)))
colData <- read.table(file = "revised_drug_response.txt", header = TRUE, sep = '\t',stringsAsFactors = F)

gene_ref <- read.csv("D:/work/tnbc_cell/mapped/Ensemble_genes_99_Human_genes(GRCh38.p13)mart_export.txt",header = T,stringsAsFactors = F, sep = "\t")
gene_ref <- gene_ref[,4:5]
gene_ref <- unique(gene_ref)

library(DESeq2); library(dplyr); library(tidyverse); library(ggplot2); library(ComplexHeatmap); library(pheatmap); library(pathfindR); library(e1071);library(org.Hs.eg.db)

# Starting Iterative steps with first drug with the
# column indice of 2 in colData 
j <- 2
repeat {
i <- names(colData[j])

  colData_tmp <- as.data.frame(colData %>%  select("Cells", all_of(i)))
  colData_tmp[,2] <- as.factor(colData_tmp[,2])
  dataset <- DESeqDataSetFromMatrix(countData = countData,colData = colData_tmp,design = formula(paste("~", i)))
  dds <- DESeq(dataset,minReplicatesForReplace = Inf)
  vsd <- vst(dds,blind = F)
  results_tmp <- results(dds, contrast=c(i, "R", "S"),cooksCutoff = T, independentFiltering=T, pAdjustMethod = "bonferroni", alpha = 0.05)
  results_tmp <- na.omit(results_tmp[order(results_tmp$padj),])
  ##### padj part #####
  # This part is for strictly adjusted p-values due to multiple testing
  results_tmp_padj <- results_tmp[results_tmp$padj <= 0.1,]
  
  # Replace Ensembl Gene Symbol IDs (ENSG000...s) to Offical Gene Symbols
  mat  <- assay(vsd[rownames(results_tmp_padj),])
  tmp <- rownames_to_column(as.data.frame(mat))
  tmp <- left_join(tmp, gene_ref,by = c( "rowname" = "Gene.stable.ID"))
  rownames(tmp) <- make.unique(tmp$Gene.name)
  mat<- tmp[,2:18]
  anno <-as.data.frame(colData(vsd)[, c("Cells", i)])
  pheatmap(mat = mat, annotation_col = anno, clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan", clustering_method = "ward.D2")
  setwd("D:/work/tnbc_cell/mapped/svm/padj")
  write.table(tmp, sep = "\t", row.names = F,col.names = T,quote = F,paste0(i, "_padj.txt"))
  setwd("D:/work/tnbc_cell/mapped/svm")
  ##### This part is for combined filtering of lfc and p-value####
  # and then applying SVM-RFE algorithm with sigFeature package
  results_tmp_p_lfc <- results_tmp[results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange <= -1.5 |results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange >= 1.5,]
   
####
  
  ##### selection with logfc and p-value filtering #####
  results_tmp_p_lfc <- na.omit(results_tmp_p_lfc[order(results_tmp_p_lfc$pvalue),])
  results_tmp_p_lfc <- results_tmp_p_lfc[1:100,]
  mat  <- assay(vsd[rownames(results_tmp_p_lfc),])
  tmp <- rownames_to_column(as.data.frame(mat))
  tmp <- left_join(tmp, gene_ref,by = c( "rowname" = "Gene.stable.ID"))
  rownames(tmp) <- make.unique(tmp$Gene.name)
  mat<- tmp[,2:18]
  anno <-as.data.frame(colData(vsd)[, c("Cells", i)])
  pheatmap(mat = mat, annotation_col = anno,cutree_rows = 4, cutree_cols = 4,fontsize_row = 7.5, clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan", clustering_method = "ward.D2")
  setwd("D:/work/tnbc_cell/mapped/svm/p_lfc")
  write.table(tmp, sep = "\t", row.names = F,col.names = T,quote = F,paste0(i, "_p_lfc.txt"))
  setwd("D:/work/tnbc_cell/mapped/svm")
  ##### pathfindR  Enrichment Part ####
  results_tmp_p_lfc <- results_tmp[results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange <= -1.5 |results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange >= 1.5,]
  tmp <- rownames_to_column(as.data.frame(results_tmp_p_lfc))
  tmp <- left_join(tmp, gene_ref,by = c( "rowname" = "Gene.stable.ID"))
  tmp <- tmp %>% select(Gene.name, log2FoldChange,pvalue)
  
  # assign(paste0("pf_",i),run_pathfindR(tmp,p_val_threshold = 0.05,visualize_enriched_terms = F,output_dir = i))

  term_gene_graph(get(paste0("pf_",i)),use_description = T)
  j <- as.numeric(match(i,names(colData)))
  j <- j+1
  if(j > 54) {
    break
  }
  
}

### This part is added to results later on ####
setwd("D:/work/tnbc_cell/mapped/svm/term_graphs/")
j <- 2
repeat {
i <- names(colData[j])
svg(paste0(i,"_term_graph.svg"),width = 8, height = 10)
term_gene_graph(get(paste0("pf_",i)),use_description = T)
dev.off()
j <- as.numeric(match(i,names(colData)))
j <- j+1
if(j >= 54) {
  break
}

}
#### ######
#### 
#### 
#### 
j <- 2
##### Run this part until you obtain the J = 54. Somehow Repeat function makes failure on svg()
##### No time to look further in this loop issue.
i <- names(colData[j])
pdf(paste0(i,"_full_term_graph.pdf"),width = 12, height = 12,family = "serif")
term_gene_graph(get(paste0("pf_",i)),use_description = T,num_terms = Inf)
dev.off()
j <- as.numeric(match(i,names(colData)))
j <- j+1



#### Lastly comparison of the 2 most KEGG pathway enriched 
#### drugs that are targetting same pathways
pdf.options(family = "serif")
# Modify "combined_results_graph" function to only highlight commonly upregulated or downregulated DEGs
########
combined_results_graph <- function (combined_df, selected_terms = "common", use_description = FALSE, 
                                    layout = "stress", node_size = "num_genes") 
{
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }
  ID_column <- ifelse(use_description, "Term_Description", 
                      "ID")
  val_node_size <- c("num_genes", "p_val")
  if (!node_size %in% val_node_size) {
    stop("`node_size` should be one of ", paste(dQuote(val_node_size), 
                                                collapse = ", "))
  }
  if (!is.data.frame(combined_df)) 
    stop("`combined_df` should be a data frame")
  necessary_cols <- c(ID_column, "combined_p", "Up_regulated_A", 
                      "Down_regulated_A", "Up_regulated_B", "Down_regulated_B")
  if (!all(necessary_cols %in% colnames(combined_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "), 
                 "must be present in `results_df`!"), collapse = " "))
  }
  if (any(selected_terms == "common")) {
    combined_df <- combined_df[combined_df$status == "common", 
                               ]
  }
  else {
    if (!any(selected_terms %in% combined_df[, ID_column])) 
      stop("None of the `selected_terms` are in the combined results!")
    combined_df <- combined_df[combined_df[, ID_column] %in% 
                                 selected_terms, ]
  }
  graph_df <- data.frame()
  for (i in base::seq_len(nrow(combined_df))) {
    up_genes <- c(unlist(strsplit(combined_df$Up_regulated_A[i], 
                                  ", ")), unlist(strsplit(combined_df$Up_regulated_B[i], 
                                                          ", ")))
    down_genes <- c(unlist(strsplit(combined_df$Down_regulated_A[i], 
                                    ", ")), unlist(strsplit(combined_df$Down_regulated_B[i], 
                                                            ", ")))
    genes <- c(up_genes, down_genes)
    genes <- genes[!is.na(genes)]
    for (gene in genes) {
      graph_df <- rbind(graph_df, data.frame(Term = combined_df[i, 
                                                                ID_column], Gene = gene, stringsAsFactors = FALSE))
    }
  }
  graph_df <- unique(graph_df)
  up_genes_A <- unlist(lapply(combined_df$Up_regulated_A, function(x) unlist(strsplit(x, 
                                                                                      ", "))))
  down_genes_A <- unlist(lapply(combined_df$Down_regulated_A, 
                                function(x) unlist(strsplit(x, ", "))))
  up_genes_B <- unlist(lapply(combined_df$Up_regulated_B, function(x) unlist(strsplit(x, 
                                                                                      ", "))))
  down_genes_B <- unlist(lapply(combined_df$Down_regulated_B, 
                                function(x) unlist(strsplit(x, ", "))))
  terms_A <- combined_df[!is.na(combined_df$lowest_p_A) & is.na(combined_df$lowest_p_B), 
                         ID_column]
  terms_B <- combined_df[is.na(combined_df$lowest_p_A) & !is.na(combined_df$lowest_p_B), 
                         ID_column]
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  igraph::V(g)$type <- ifelse(names(igraph::V(g)) %in% terms_A, 
                              "A-only term", ifelse(names(igraph::V(g)) %in% 
                                                      terms_B, "B-only term", ifelse(names(igraph::V(g)) %in% 
                                                                                       combined_df[, ID_column], "common term", "gene")))
  if (node_size == "num_genes") {
    sizes <- igraph::degree(g)
    sizes <- ifelse(grepl("term", igraph::V(g)$type), 
                    sizes, 2)
    size_label <- "# genes"
  }
  else {
    idx <- match(names(igraph::V(g)), combined_df[, ID_column])
    sizes <- -log10(combined_df$combined_p[idx])
    sizes[is.na(sizes)] <- 2
    size_label <- "-log10(p)"
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.5
  igraph::V(g)$frame.color <- "gray"
  cond_up_A <- names(igraph::V(g)) %in% up_genes_A
  cond_up_B <- names(igraph::V(g)) %in% up_genes_B
  cond_down_A <- names(igraph::V(g)) %in% down_genes_A
  cond_down_B <- names(igraph::V(g)) %in% down_genes_B
  missing_A <- !cond_up_A & !cond_down_A
  missing_B <- !cond_up_B & !cond_down_B
  up_cond <- (cond_up_A & cond_up_B) #| (missing_A & cond_up_B) | 
  #(cond_up_A & missing_B)
  down_cond <- (cond_down_A & cond_down_B) #| (missing_A & cond_down_B) | 
  #(cond_down_A & missing_B)
  igraph::V(g)$for_coloring <- ifelse(igraph::V(g)$type == 
                                        "common term", "Common term", ifelse(igraph::V(g)$type == 
                                                                               "A-only term", "A-only term", ifelse(igraph::V(g)$type == 
                                                                                                                      "B-only term", "B-only term", ifelse(up_cond, 
                                                                                                                                                           "Up gene", ifelse(down_cond, "Down gene", 
                                                                                                                                                                             "Conflict or NPiB")))))
  p <- ggraph::ggraph(g, layout = layout)
  p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey")
  p <- p + ggraph::geom_node_point(ggplot2::aes_(color = ~for_coloring, 
                                                 size = ~size))
  p <- p + ggplot2::scale_size(range = c(5, 10), breaks = round(seq(round(min(igraph::V(g)$size)), 
                                                                    round(max(igraph::V(g)$size)), length.out = 4)), name = size_label)
  p <- p + ggplot2::theme_void()
  p <- p + ggraph::geom_node_text(ggplot2::aes_(label = ~name), 
                                  nudge_y = 0.2)
  vertex_cols <- c(`Common term` = "#FCCA46", `A-only term` = "#9FB8AD", 
                   `B-only term` = "#619B8A", `Up gene` = "green", 
                   `Down gene` = "red", `Conflict or NPiB` = "gray")
  p <- p + ggplot2::scale_colour_manual(values = vertex_cols, 
                                        name = NULL)
  # p <- p + ggplot2::ggtitle("Combined Terms Graph")
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(p)
}


##############

combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_ABT737,result_B = pf_Navitoclax),use_description = T) # Apoptosis Regulation
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_RO.3306,result_B = pf_Wee1.Inhibitor),use_description = T) # Cell Cycle
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Entinostat,result_B = pf_Vorinostat),use_description = T) # Chromatin histone acetylation
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_BDP.00009066,result_B = pf_GSK269962A),use_description = T) # Cytoskeleton
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Epirubicin,result_B = pf_Irinotecan),use_description = T) # DNA Replication
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Alisertib,result_B = pf_Tozasertib),use_description = T) # Mitosis
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Sapitinib,result_B = pf_Afatinib),use_description = T) # EGFR Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Afuresertib,result_B = pf_Ipatasertib),use_description = T) # PI3K/MTOR Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Trametinib,result_B = pf_VX.11e),use_description = T) # RTK Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_NU7441,result_B = pf_VE821),use_description = T) # Genome_Integrity
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_BMS.536924,result_B = pf_NVP.ADW742),use_description = T) # IGF1R Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Dactolisib,result_B = pf_AZD8186),use_description = T)
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Crizotinib,result_B = pf_AZD4547),use_description = T)
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_AZ6102,result_B = pf_WIKI4),use_description = T)

##  intersection of the most similar drugs
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_AZD4547,result_B = pf_ZM447439),use_description = T)
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_AZD6738,result_B = pf_Pevonedistat),use_description = T)
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_NVP.ADW742,result_B = pf_WIKI4),use_description = T)
