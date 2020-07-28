setwd("D:/work/tnbc_cell/mapped/svm")
countData <- read.table(file = "whole_exp.txt" , header = TRUE, sep = '\t',row.names = 1)
countData <- countData[,sort(names(countData))]
names(countData) <-(gsub("\\.", "\\-", names(countData)))
colData <- read.table(file = "revised_drug_response.txt", header = TRUE, sep = '\t',stringsAsFactors = F)

gene_ref <- read.csv("D:/work/tnbc_cell/mapped/Ensemble_genes_99_Human_genes(GRCh38.p13)mart_export.txt",header = T,stringsAsFactors = F, sep = "\t")
gene_ref <- gene_ref[,4:5]
gene_ref <- unique(gene_ref)

library(DESeq2); library(dplyr); library(tidyverse); library(ggplot2); library(ComplexHeatmap); library(pheatmap); library(pathfindR); library(e1071);library(org.Hs.eg.db); #library(sigFeature)

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
  # # 
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
# library(svglite)
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

combined_results_graph(node_size = "p_val",use_description = T,combine_pathfindR_results(result_A = pf_ABT737,result_B = pf_Navitoclax)) # Apoptosis Regulation


combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_RO.3306,result_B = pf_Wee1.Inhibitor),use_description = T) # Cell Cycle
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Entinostat,result_B = pf_Vorinostat),use_description = T) # Chromatin histone acetylation
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_BDP.00009066,result_B = pf_GSK269962A),use_description = T) # Cytoskeleton
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Epirubicin,result_B = pf_Irinotecan),use_description = T) # DNA Replication
# combined_results_graph(combine_pathfindR_results(result_A = pf_Epirubicin,result_B = pf_Gemcitabine),use_description = T) # DNA Replication
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Alisertib,result_B = pf_Tozasertib),use_description = T) # Mitosis
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Sapitinib,result_B = pf_Afatinib),use_description = T) # EGFR Signalling
# combined_results_graph(combine_pathfindR_results(result_A = pf_WIKI4,result_B = pf_AZ6102),use_description = T) # WNT Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Afuresertib,result_B = pf_Ipatasertib),use_description = T) # PI3K/MTOR Signalling
combined_results_graph(combine_pathfindR_results(result_A = pf_Vinorelbine,result_B = pf_Vinblastine),use_description = T) # Mitosis
combined_results_graph(combine_pathfindR_results(result_A = pf_Cediranib,result_B = pf_Crizotinib),use_description = T) # RTK Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Trametinib,result_B = pf_VX.11e),use_description = T) # RTK Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_NU7441,result_B = pf_VE821),use_description = T) # Genome_Integrity
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_BMS.536924,result_B = pf_NVP.ADW742),use_description = T) # IGF1R Signalling
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Dactolisib,result_B = pf_AZD8186),use_description = T) 
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_Crizotinib,result_B = pf_AZD4547),use_description = T)
combined_results_graph(node_size = "p_val",combine_pathfindR_results(result_A = pf_AZ6102,result_B = pf_WIKI4),use_description = T)




