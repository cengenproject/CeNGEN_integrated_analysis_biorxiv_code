library(circlize)
library(ComplexHeatmap)


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

consensus_Recall_matrix <- data.frame(matrix(data = 0,
                                                  nrow = length(unique(temp_mat[,1])),
                                                  ncol = length(unique(temp_mat[,1])),
                                                  dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(auc_TPR_consensus_p.hmp_directional)){
    consensus_Recall_matrix[cell1, cell2] <- auc_TPR_consensus_p.hmp_directional[[contrast]]
  } else {consensus_Recall_matrix[cell1, cell2] <- 0} }

for(cell in rownames(consensus_Recall_matrix)){
  consensus_Recall_matrix[cell, cell] <- NA
}

raw_Recall_matrix <- data.frame(matrix(data = 0,
                                            nrow = length(unique(temp_mat[,1])),
                                            ncol = length(unique(temp_mat[,1])),
                                            dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(auc_TPR_raw_PValue_directional)){
    raw_Recall_matrix[cell1, cell2] <- auc_TPR_raw_PValue_directional[[contrast]]
  } else {raw_Recall_matrix[cell1, cell2] <- 0}}

for(cell in rownames(raw_Recall_matrix)){
  raw_Recall_matrix[cell, cell] <- NA
}

consensus_Recall_matrix <- consensus_Recall_matrix[rowSums(consensus_Recall_matrix, na.rm = T) > 0, colSums(consensus_Recall_matrix, na.rm = T) >0]
raw_Recall_matrix <- raw_Recall_matrix[rowSums(raw_Recall_matrix, na.rm = T) > 0, colSums(raw_Recall_matrix, na.rm = T) >0]
dim(raw_Recall_matrix)
dim(consensus_Recall_matrix)

col_fun <- circlize::colorRamp2(breaks = c(0, max(unlist(raw_Recall_matrix), na.rm = T)), c('#DEDCDC', '#AF1302'))

Recall_diff_mat <- consensus_Recall_matrix - raw_Recall_matrix

##Figure SCi
pdf("raw_Recall_matrix.pdf",width=5.5,height=5)
Heatmap(raw_Recall_matrix,  cluster_rows = F, cluster_columns = F,
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()  

##Figure SCii
pdf("consensus_Recall_matrix.pdf",width=5.5,height=5)
Heatmap(consensus_Recall_matrix,  cluster_rows = F, cluster_columns = F, 
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()  

col_fun2 <- circlize::colorRamp2(breaks = c(-1, 0, 1), c('#1e3dab', '#F7f7f7', '#AF1302'))

##Figure SCiii
pdf("Recall_diff_mat.pdf",width=5.5,height=5)
Heatmap(Recall_diff_mat,  cluster_rows = F, cluster_columns = F,
        col = col_fun2, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()  


Recall_for_histogram.df <- data.frame(raw_Recall = unlist(raw_Recall_matrix),
                                              consensus_Recall =unlist(consensus_Recall_matrix),
                                              Difference_Recall = unlist(Recall_diff_mat))
Recall_for_histogram.df <- na.omit(Recall_for_histogram.df)

##Figure SCiv
ggplot(Recall_for_histogram.df) + 
  stat_density(aes(x = Difference_Recall, y = ..count..), fill = '#8940c6', color = 'black', alpha = 0.6, adjust = 1) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) +
  xlim(-1, 1) +
  xlab('Difference in Differential Expression Recall\nIntegrated minus Bulk')
ggsave('Recall_for_histogram.df_density.pdf', width = 6, height = 6)
