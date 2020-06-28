setwd("D:/work/tnbc_cell/mapped/svm")
countData <- read.table(file = "whole_exp.txt" , header = TRUE, sep = '\t',row.names = 1)
countData <- countData[,sort(names(countData))]
names(countData) <-(gsub("\\.", "\\-", names(countData)))
colData <- read.table(file = "revised_drug_response.txt", header = TRUE, sep = '\t',stringsAsFactors = F)

gene_ref <- read.csv("D:/work/tnbc_cell/mapped/Ensemble_genes_99_Human_genes(GRCh38.p13)mart_export.txt",header = T,stringsAsFactors = F, sep = "\t")
gene_ref <- gene_ref[,4:5]
gene_ref <- unique(gene_ref)

library(DESeq2); library(dplyr); library(tidyverse); library(ggplot2); library(ComplexHeatmap); library(pheatmap); library(pathfindR); library(sigFeature); library(e1071);library(org.Hs.eg.db);

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
  pheatmap(mat = mat, annotation_col = anno, clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan", clustering_method = "ward.D2", main = paste("DEGs (p_adj <= 0.1) in Res. vs Sens. cells for",i))
  setwd("D:/work/tnbc_cell/mapped/svm/padj")
  write.table(tmp, sep = "\t", row.names = F,col.names = T,quote = F,paste0(i, "_padj.txt"))
  setwd("D:/work/tnbc_cell/mapped/svm")
  ##### This part is for combined filtering of lfc and p-value####
  # and then applying SVM-RFE algorithm with sigFeature package
  results_tmp_p_lfc <- results_tmp[results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange <= -1.5 |results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange >= 1.5,]
    #### sigFeature part####
  # x <- as.data.frame(t(countData[rownames(countData) %in% rownames(results_tmp_p_lfc),]))
  # temp <- rownames_to_column(as.data.frame(t(x)))
  # temp <- left_join(temp, gene_ref,by = c( "rowname" = "Gene.stable.ID"))
  # rownames(temp) <- make.unique(temp$Gene.name)
  # x <- as.data.frame(t(temp[,2:18]))
  # y <- as.factor(colData_tmp[,2])
  # 
  # pvals <- sigFeaturePvalue(x,y)
  # 
  # system.time(sigfeatureRankedList <- sigFeature(x, y))
  # 
  # pheatmap(x[ ,sigfeatureRankedList[1:100]], scale="row", 
  #        clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan", clustering_method = "ward.D2",main = paste("sigFeature results Res. vs Sens.for",i))
  # 
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
  pheatmap(mat = mat, annotation_col = anno,cutree_rows = 4, cutree_cols = 4,fontsize_row = 7.5, clustering_distance_rows = "manhattan",clustering_distance_cols = "manhattan", clustering_method = "ward.D2", main = paste("Top 100 DEGs (p-value <= 0.05 & lfc +- 1.5) in Res. vs Sens. cells for",i))
  setwd("D:/work/tnbc_cell/mapped/svm/p_lfc")
  write.table(tmp, sep = "\t", row.names = F,col.names = T,quote = F,paste0(i, "_p_lfc.txt"))
  setwd("D:/work/tnbc_cell/mapped/svm")
  ##### pathfindR  Enrichment Part ####
  results_tmp_p_lfc <- results_tmp[results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange <= -1.5 |results_tmp$pvalue <= 0.05 & results_tmp$log2FoldChange >= 1.5,]
  tmp <- rownames_to_column(as.data.frame(results_tmp_p_lfc))
  tmp <- left_join(tmp, gene_ref,by = c( "rowname" = "Gene.stable.ID"))
  tmp <- tmp %>% select(Gene.name, log2FoldChange,pvalue)
  
  assign(paste0("pf_",i),run_pathfindR(tmp,p_val_threshold = 0.05,visualize_enriched_terms = F,output_dir = i))
  
  cluster_enriched_terms(get(paste0("pf_",i)),use_description = T)
  
  j <- as.numeric(match(i,names(colData)))
  j <- j+1
  if(j > 53) {
    break
  }
  
}
#### Can be modified for later usage ####
# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from=plots.png.paths, to="D:/work/tnbc_cell/mapped/svm/plots")


