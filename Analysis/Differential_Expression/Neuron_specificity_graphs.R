


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

consensus_Specificity_matrix <- data.frame(matrix(data = 0,
                                         nrow = length(unique(temp_mat[,1])),
                                         ncol = length(unique(temp_mat[,1])),
                                         dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(auc_FPR_consensus_p.hmp_directional)){
    consensus_Specificity_matrix[cell1, cell2] <- 1-auc_FPR_consensus_p.hmp_directional[[contrast]]
  } else {consensus_Specificity_matrix[cell1, cell2] <- 0} }

for(cell in rownames(consensus_Specificity_matrix)){
  consensus_Specificity_matrix[cell, cell] <- NA
}


harmonic_Specificity_matrix <- data.frame(matrix(data = 0,
                                        nrow = length(unique(temp_mat[,1])),
                                        ncol = length(unique(temp_mat[,1])),
                                        dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(auc_FPR_harmonized_p.hmp_directional)){
    harmonic_Specificity_matrix[cell1, cell2] <- 1-auc_FPR_harmonized_p.hmp_directional[[contrast]]
  } else {harmonic_Specificity_matrix[cell1, cell2] <- 0} }

for(cell in rownames(harmonic_Specificity_matrix)){
  harmonic_Specificity_matrix[cell, cell] <- NA
}


raw_Specificity_matrix <- data.frame(matrix(data = 0,
                                   nrow = length(unique(temp_mat[,1])),
                                   ncol = length(unique(temp_mat[,1])),
                                   dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(auc_FPR_raw_PValue_directional)){
    raw_Specificity_matrix[cell1, cell2] <- 1-auc_FPR_raw_PValue_directional[[contrast]]
  } else {raw_Specificity_matrix[cell1, cell2] <- 0}}

for(cell in rownames(raw_Specificity_matrix)){
  raw_Specificity_matrix[cell, cell] <- NA
}

harmonic_Specificity_matrix <- harmonic_Specificity_matrix[rowSums(harmonic_Specificity_matrix, na.rm = T) > 0, colSums(harmonic_Specificity_matrix, na.rm = T) >0]
consensus_Specificity_matrix <- consensus_Specificity_matrix[rowSums(consensus_Specificity_matrix, na.rm = T) > 0, colSums(consensus_Specificity_matrix, na.rm = T) >0]
raw_Specificity_matrix <- raw_Specificity_matrix[rowSums(raw_Specificity_matrix, na.rm = T) > 0, colSums(raw_Specificity_matrix, na.rm = T) >0]
dim(harmonic_Specificity_matrix)
dim(raw_Specificity_matrix)
dim(consensus_Specificity_matrix)
library(circlize)

col_fun <- circlize::colorRamp2(breaks = c(min(unlist(raw_Specificity_matrix), na.rm = T), 1), c('white', '#AF1302'))
#col_fun <- circlize::colorRamp2(breaks = c(0, 1), c('white', '#AF1302'))


Specificity_diff_mat <- consensus_Specificity_matrix - raw_Specificity_matrix


### Figure SDi
pdf("raw_Specificity_matrix.pdf",width=5.5,height=5)
Heatmap(raw_Specificity_matrix,  cluster_rows = F, cluster_columns = F,
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()

### Figure SDii
pdf("consensus_Specificity_matrix.pdf",width=5.5,height=5)
Heatmap(consensus_Specificity_matrix,  cluster_rows = F, cluster_columns = F, 
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()

col_fun2 <- circlize::colorRamp2(breaks = c(-0.5, 0, 0.5), c('#1e3dab', '#F7f7f7', '#AF1302'))


### Figure SDiii
pdf("Specificity_diff_mat.pdf",width=5.5,height=5)
Heatmap(Specificity_diff_mat,  cluster_rows = F, cluster_columns = F,
        col = col_fun2, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()



Specificity_for_histogram.df <- data.frame(raw_Specificity = unlist(raw_Specificity_matrix),
                                              consensus_Specificity =unlist(consensus_Specificity_matrix),
                                              Difference_in_Specificity = unlist(Specificity_diff_mat))
Specificity_for_histogram.df <- na.omit(Specificity_for_histogram.df)


### Figure SDiv
ggplot(Specificity_for_histogram.df) + 
  geom_density(aes(x = Difference_in_Specificity), fill = '#8940c6', color = 'black', alpha = 0.6) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) +
  xlim(-0.35, 0.35) +
  xlab('Difference in Differential Expression Specificity\nIntegrated minus Bulk')
ggsave('Specificity_for_histogram.df_density.pdf', width = 6, height = 6)



