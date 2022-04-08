### this file builds directly from the outputs of "Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R" and requires running the differential expression analysis and that script.


### build an empty matrix of the correct size
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



matched_contrast

consensus_Accuracy_matrix <- data.frame(matrix(data = 0,
                                             nrow = length(unique(temp_mat[,1])),
                                             ncol = length(unique(temp_mat[,1])),
                                             dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(consensus_p.hmp_accuracy)){
    consensus_Accuracy_matrix[cell1, cell2] <- consensus_p.hmp_accuracy[[contrast]]
  } else {consensus_Accuracy_matrix[cell1, cell2] <- 0} }

for(cell in rownames(consensus_Accuracy_matrix)){
  consensus_Accuracy_matrix[cell, cell] <- NA
}


raw_Accuracy_matrix <- data.frame(matrix(data = 0,
                                       nrow = length(unique(temp_mat[,1])),
                                       ncol = length(unique(temp_mat[,1])),
                                       dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(raw_PValue_accuracy)){
    raw_Accuracy_matrix[cell1, cell2] <- raw_PValue_accuracy[[contrast]]
  } else {raw_Accuracy_matrix[cell1, cell2] <- 0}}

for(cell in rownames(raw_Accuracy_matrix)){
  raw_Accuracy_matrix[cell, cell] <- NA
}

consensus_Accuracy_matrix <- consensus_Accuracy_matrix[rowSums(consensus_Accuracy_matrix, na.rm = T) > 0, colSums(consensus_Accuracy_matrix, na.rm = T) >0]
raw_Accuracy_matrix <- raw_Accuracy_matrix[rowSums(raw_Accuracy_matrix, na.rm = T) > 0, colSums(raw_Accuracy_matrix, na.rm = T) >0]
dim(harmonic_Accuracy_matrix)
dim(raw_Accuracy_matrix)
dim(consensus_Accuracy_matrix)
library(circlize)

col_fun <- circlize::colorRamp2(breaks = c(0, 0.7, max(unlist(raw_Accuracy_matrix), na.rm = T)), c('white','lightgrey', '#AF1302'))

col_fun <- circlize::colorRamp2(breaks = c(min(raw_Accuracy_matrix, na.rm = T), 1), c('white', '#AF1302'))


Accuracy_diff_mat <- consensus_Accuracy_matrix - raw_Accuracy_matrix

row_dend = dendsort(hclust(dist((Accuracy_diff_mat))))

### Figure S3Ei
pdf("raw_Accuracy_matrix.pdf",width=5.5,height=5)
Heatmap(raw_Accuracy_matrix,  cluster_rows = F, cluster_columns = F,
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()  

### Figure S3Eii
pdf("consensus_Accuracy_matrix.pdf",width=5.5,height=5)
Heatmap(consensus_Accuracy_matrix,  cluster_rows = F, cluster_columns = F, 
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()  

col_fun2 <- circlize::colorRamp2(breaks = c(-0.5, 0, 0.5), c('#1e3dab', '#F7f7f7', '#AF1302'))

### Figure S3Eiii
pdf("Accuracy_diff_mat.pdf",width=5.5,height=5)
Heatmap(Accuracy_diff_mat,  cluster_rows = F, cluster_columns = F,
        col = col_fun2, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()  



quantile(unlist(Accuracy_diff_mat), seq(0,1,0.01), na.rm = T) %>% data.frame(difference = .,
                                                                           percents = seq(0,1,0.01)) %>%
  ggplot() + geom_hline(yintercept = 0) + geom_step(aes(x = percents, y = difference), color = '#AF1302', size = 1) +
  theme_classic(base_size = 35) +
  xlab('Percentile: 0 -> 100') + ylab('Difference in Accuracy\nIntegrated minus Raw') +
  ylim(-0.5,0.5)



Accuracy_for_histogram.df <- data.frame(raw_accuracy = unlist(raw_Accuracy_matrix),
                                  consensus_accuracy =unlist(consensus_Accuracy_matrix),
                                  difference_in_accuracy = unlist(Accuracy_diff_mat))
Accuracy_for_histogram.df <- na.omit(Accuracy_for_histogram.df)


### Figure 3a
ggplot(Accuracy_for_histogram.df) + 
  geom_density(aes(x = raw_accuracy), fill = '#4ebfde', color = 'black', alpha = 0.6) +
  geom_density(aes(x = consensus_accuracy), fill = '#C8594b', color = 'black', alpha = 0.6) +
  theme_classic(base_size = 20) +
  xlim(0.5,1) + 
  xlab('Pairwise Differential Expression Accuracy')
ggsave('Accuracy_for_histogram.df_density_both_113021.pdf', units = 'in', width = 10, height = 6.55)

### Figure 3b
ggplot(Accuracy_for_histogram.df) + 
  geom_density(aes(x = difference_in_accuracy), fill = '#8940c6', color = 'black', alpha = 0.6) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) +
  xlim(-0.4,0.4) +
  xlab('Difference in Differenti')
ggsave('Accuracy_for_histogram.df_density_difference_113021.pdf', units = 'in', width = 10, height = 6.55)



