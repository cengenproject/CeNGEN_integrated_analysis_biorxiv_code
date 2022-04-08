### this file builds directly from the outputs of _____ and requires running the differential expression analysis and that script.



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

unique(temp_mat[,1])


matched_contrast


## duplicate the matrix and fill in the appropriate values from the MCC list ----
consensus_MCC_matrix <- data.frame(matrix(data = 0,
                                                  nrow = length(unique(temp_mat[,1])),
                                                  ncol = length(unique(temp_mat[,1])),
                                                  dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(consensus_MCC)){
    consensus_MCC_matrix[cell1, cell2] <- consensus_MCC[[contrast]]
  } else {consensus_MCC_matrix[cell1, cell2] <- 0} }



for(cell in rownames(consensus_MCC_matrix)){
  consensus_MCC_matrix[cell, cell] <- NA
}


## duplicate the matrix and fill in the appropriate values from the bulk MCC list ----

raw_MCC_matrix <- data.frame(matrix(data = 0,
                                            nrow = length(unique(temp_mat[,1])),
                                            ncol = length(unique(temp_mat[,1])),
                                            dimnames = list(unique(temp_mat[,1]), unique(temp_mat[,1]))))


for(contrast in rownames(temp_mat)){
  cell1 <- temp_mat[contrast,1]
  cell2 <- temp_mat[contrast,2]
  
  if(contrast %in% names(raw_pValue_MCC)){
    raw_MCC_matrix[cell1, cell2] <- raw_pValue_MCC[[contrast]]
  } else {raw_MCC_matrix[cell1, cell2] <- 0}}

for(cell in rownames(raw_MCC_matrix)){
  raw_MCC_matrix[cell, cell] <- NA
}

consensus_MCC_matrix <- consensus_MCC_matrix[rowSums(consensus_MCC_matrix, na.rm = T) > 0, colSums(consensus_MCC_matrix, na.rm = T) >0]
raw_MCC_matrix <- raw_MCC_matrix[rowSums(raw_MCC_matrix, na.rm = T) > 0, colSums(raw_MCC_matrix, na.rm = T) >0]
dim(raw_MCC_matrix)
dim(consensus_MCC_matrix)
library(circlize)

col_fun <- circlize::colorRamp2(breaks = c(min(unlist(raw_MCC_matrix), na.rm = T), 1), c('white', '#AF1302'))
col_fun <- circlize::colorRamp2(breaks = c(-1, 0, 1), c('blue', 'white', '#AF1302'))


MCC_diff_mat <- consensus_MCC_matrix - raw_MCC_matrix

row_dend = dendsort(hclust(dist((Recall_diff_mat))))


##Figure SFi
pdf("raw_MCC_matrix.pdf",width=5.5,height=5)
Heatmap(raw_MCC_matrix,  cluster_rows = F, cluster_columns = F,
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()

##Figure SFii
pdf("consensus_MCC_matrix.pdf",width=5.5,height=5)
Heatmap(consensus_MCC_matrix,  cluster_rows = F, cluster_columns = F, 
        col = col_fun, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()


quantile(unlist(MCC_diff_mat), na.rm = T)
col_fun2 <- circlize::colorRamp2(breaks = c(-0.85, 0, 0.85), c('#1e3dab', '#F7f7f7', '#AF1302'))

### Figure S3Fiii
pdf("MCC_diff_mat.pdf",width=5.5,height=5)
Heatmap(MCC_diff_mat,  cluster_rows = F, cluster_columns = F,
        col = col_fun2, show_row_dend = F, show_column_dend = F,
        name = ' ', row_names_side = 'left', use_raster = T)
dev.off()

quantile(unlist(MCC_diff_mat), seq(0,1,0.01), na.rm = T)

quantile(unlist(MCC_diff_mat), seq(0,1,0.01), na.rm = T) %>% data.frame(difference = .,
                                                                                percents = seq(0,1,0.01)) %>%
  ggplot() + geom_hline(yintercept = 0) + geom_step(aes(x = percents, y = difference), color = '#AF1302', size = 3) +
  theme_classic(base_size = 20) +
  xlab('Percentile: 0 -> 100') + ylab('Difference in MCC (True Negative Rate)\nIntegrated minus Raw') +
  ylim(-0.5,0.5)


MCC_for_histogram.df <- data.frame(raw_MCC = unlist(raw_MCC_matrix),
                                           consensus_MCC =unlist(consensus_MCC_matrix),
                                           Difference_in_MCC = unlist(MCC_diff_mat))
MCC_for_histogram.df <- na.omit(MCC_for_histogram.df)



### Figure 3c
ggplot(MCC_for_histogram.df) + 
  geom_density(aes(x = raw_MCC), fill = '#4ebfde', color = 'black', alpha = 0.6) +
  geom_density(aes(x = consensus_MCC), fill = '#C8594b', color = 'black', alpha = 0.6) +
  theme_classic(base_size = 20) +
  xlim(-0.5,1) + 
  xlab('Bulk Differential Expression MCC')
ggsave('MCC_for_histogram.df_density.pdf', units = 'in', width = 10, height = 6.55)



### spiky histogram version
ggplot(MCC_for_histogram.df) + 
  geom_area(aes(x = Difference_in_MCC), fill = '#8940c6', color = 'black', alpha = 0.6,  stat = "bin", bins = 30) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) +
  xlim(-1, 1) +
  xlab('DEX MCC (Int minus Bulk')

### Figure 3d & S3Fiv
ggplot(MCC_for_histogram.df) + 
  geom_density(aes(x = Difference_in_MCC), fill = '#8940c6', color = 'black', alpha = 0.6) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
  theme_classic(base_size = 20) +
  xlim(-1, 1) +
  xlab('DEX MCC (Int minus Bulk')
ggsave('MCC_for_histogram._difference.df_density.pdf', width = 10, height = 6.55)
ggsave('MCC_for_histogram._difference.df_density_small.pdf', width = 5, height = 5)

