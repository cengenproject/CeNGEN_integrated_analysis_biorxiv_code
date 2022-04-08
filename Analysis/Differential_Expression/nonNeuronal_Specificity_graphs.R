library(circlize)
library(ComplexHeatmap)



non_neuronal_gt_slim <- rep(0, length(nn_genes))
names(non_neuronal_gt_slim) <- nn_genes
non_neuronal_gt_slim <- non_neuronal_gt_slim[intersect(names(non_neuronal_gt_slim), rownames(many_harmonized$`ADL-(AFD)`))]


auc_nn_FPR_consensus_p.hmp_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- many_harmonized[[contrast]]
    
    testing_gt1 <- non_neuronal_gt_slim
    
    
    testing_harmony_true <- testing_harmony[names(testing_gt1),3]
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1), 2] < 2] <- NA
    
    FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = 40)
    
    return(FPR)
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- many_harmonized[[match]]
      
      testing_gt1 <- non_neuronal_gt_slim
      
      
      testing_harmony_true <- testing_harmony[names(testing_gt1),3]
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1), 2] > (-2)] <- NA
      
      
      FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = 40)
      
      return(FPR)
      
    }
    else{NA}
  }
})
auc_nn_FPR_raw_PValue_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- raw_qlfs[[contrast]]
    
    testing_gt1 <- non_neuronal_gt_slim
    
    
    testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] < 2] <- NA
    
    FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
    
    return(FPR)
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- raw_qlfs[[match]]
      
      testing_gt1 <- non_neuronal_gt_slim
      
      
      testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] > (-2)] <- NA
      
      FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
      
      return(FPR)
      
    }
    else{NA}
  }
})
auc_nn_FPR_consensus_p.hmp_directional <- na.omit(auc_nn_FPR_consensus_p.hmp_directional)
auc_nn_FPR_raw_PValue_directional <- na.omit(auc_nn_FPR_raw_PValue_directional)

summary(auc_nn_FPR_raw_PValue_directional)
summary(auc_nn_FPR_consensus_p.hmp_directional)



temp_mat <- lapply(matched_contrast[,1], function(contrasts){
  cell1 <- str_split_fixed(contrasts, '-', 2)[,1]
  cell2 <- str_split_fixed(contrasts, '-', 2)[,2]
  cell2 <- gsub(pattern = '\\(|\\)', replacement = '', x = cell2)
  
  #print(c(cell1, cell2))
  return(c(cell1, cell2))
  #return(cell2)
})

temp_mat <- do.call(rbind, temp_mat)
rownames(temp_mat) <- matched_contrast[,1]

unique(temp_mat[,1])


matched_contrast

consensus_nn_Specificity_matrix <- data.frame(matrix(data = 0,
                                                     nrow = length(unique(temp_mat[,1])),
                                                     ncol = length(unique(temp_mat[,1])),
                                                     dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(auc_nn_FPR_consensus_p.hmp_directional)){
    consensus_nn_Specificity_matrix[cell1, cell2] <- 1-auc_nn_FPR_consensus_p.hmp_directional[[contrast]]
  } else {consensus_nn_Specificity_matrix[cell1, cell2] <- 0} }

for(cell in rownames(consensus_nn_Specificity_matrix)){
  consensus_nn_Specificity_matrix[cell, cell] <- NA
}




raw_nn_Specificity_matrix <- data.frame(matrix(data = 0,
                                               nrow = length(unique(temp_mat[,1])),
                                               ncol = length(unique(temp_mat[,1])),
                                               dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(auc_nn_FPR_raw_PValue_directional)){
    raw_nn_Specificity_matrix[cell1, cell2] <- 1-auc_nn_FPR_raw_PValue_directional[[contrast]]
  } else {raw_nn_Specificity_matrix[cell1, cell2] <- 0}}

for(cell in rownames(raw_nn_Specificity_matrix)){
  raw_nn_Specificity_matrix[cell, cell] <- NA
}

consensus_nn_Specificity_matrix <- consensus_nn_Specificity_matrix[rowSums(consensus_nn_Specificity_matrix, na.rm = T) > 0, colSums(consensus_nn_Specificity_matrix, na.rm = T) >0]
raw_nn_Specificity_matrix <- raw_nn_Specificity_matrix[rowSums(raw_nn_Specificity_matrix, na.rm = T) > 0, colSums(raw_nn_Specificity_matrix, na.rm = T) >0]
dim(raw_nn_Specificity_matrix)
dim(consensus_nn_Specificity_matrix)

col_fun <- circlize::colorRamp2(breaks = c(min(unlist(raw_nn_Specificity_matrix), na.rm = T), 1), c('white', '#AF1302'))
col_fun4 <- circlize::colorRamp2(breaks = c(0.6, 1), c('white', '#AF1302'))

nn_Specificity_diff_mat <- consensus_nn_Specificity_matrix - raw_nn_Specificity_matrix

pdf("raw_nn_Specificity_matrix.pdf",width=5.5,height=5)
Heatmap(raw_nn_Specificity_matrix,  cluster_rows = F, cluster_columns = F,
        col = col_fun4, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()


pdf("consensus_nn_Specificity_matrix.pdf",width=5.5,height=5)
Heatmap(consensus_nn_Specificity_matrix,  cluster_rows = F, cluster_columns = F, 
        col = col_fun4, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()

min(raw_nn_Specificity_matrix, na.rm = T)
min(consensus_nn_Specificity_matrix, na.rm = T)


quantile(unlist(nn_Specificity_diff_mat), na.rm = T)
col_fun3 <- circlize::colorRamp2(breaks = c(-0.3, 0, 0.3), c('#1e3dab', '#F7f7f7', '#AF1302'))


pdf("nn_Specificity_diff_mat.pdf",width=5.5,height=5)
Heatmap(nn_Specificity_diff_mat,  cluster_rows = F, cluster_columns = F,
        col = col_fun3, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()



nn_Specificity_for_histogram.df <- data.frame(raw_nn_Specificity = unlist(raw_nn_Specificity_matrix),
                                              consensus_nn_Specificity =unlist(consensus_nn_Specificity_matrix),
                                              Difference_in_nn_Specificity = unlist(nn_Specificity_diff_mat))
nn_Specificity_for_histogram.df <- na.omit(nn_Specificity_for_histogram.df)



ggplot(nn_Specificity_for_histogram.df) + 
  geom_density(aes(x = Difference_in_nn_Specificity), fill = '#8940c6', color = 'black', alpha = 0.6) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) +
  xlim(-0.3, 0.3) +
  xlab('Difference in Differential Expression\nNon-neuronal SpecificityIntegrated minus Bulk')
ggsave('nn_Specificity_for_histogram.df_diff_density.pdf', width = 6, height = 6)

