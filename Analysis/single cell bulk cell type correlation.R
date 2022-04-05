### detection of single cell genes in bulk
library(dendextend)
library(stringr)
library(ComplexHeatmap)

liberal <- cengenDataSC::cengen_sc_1_bulk
medium <- cengenDataSC::cengen_sc_2_bulk
conservative <- cengenDataSC::cengen_sc_3_bulk
strict <- cengenDataSC::cengen_sc_4_bulk



genes_expressed_sc_1_list <- sapply(colnames(liberal), function(cell){
  rownames(liberal[liberal[,cell] > 0,])
})
liberal_bulk_detection <- sapply(colnames(aggr_raw_GeTMM), function(cell){
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]
  common.genes <- intersect(rownames(bulk_raw_GeTMM), genes_expressed_sc_1_list[[cell]])
  
  bulk <- bulk_raw_GeTMM[common.genes, samples_cell_type %in% c(cell)]
  
  bulk_bin <- bulk > 5
  bulk_bin_ave <- rowMeans(bulk_bin)
  return(sum(bulk_bin_ave > 0.65)/length(common.genes))
  
  
})
genes_expressed_sc_2_list <- sapply(colnames(medium), function(cell){
  rownames(medium[medium[,cell] > 0,])
})
medium_bulk_detection <- sapply(colnames(aggr_raw_GeTMM), function(cell){
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]
  common.genes <- intersect(rownames(bulk_raw_GeTMM), genes_expressed_sc_2_list[[cell]])

  bulk <- bulk_raw_GeTMM[common.genes, samples_cell_type %in% c(cell)]
  
  bulk_bin <- bulk > 5
  bulk_bin_ave <- rowMeans(bulk_bin)
  return(sum(bulk_bin_ave > 0.65)/length(common.genes))


})
genes_expressed_sc_3_list <- sapply(colnames(conservative), function(cell){
  rownames(conservative[conservative[,cell] > 0,])
})
conservative_bulk_detection <- sapply(colnames(aggr_raw_GeTMM), function(cell){
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]
  common.genes <- intersect(rownames(bulk_raw_GeTMM), genes_expressed_sc_3_list[[cell]])
  
  bulk <- bulk_raw_GeTMM[common.genes, samples_cell_type %in% c(cell)]
  
  bulk_bin <- bulk > 5
  bulk_bin_ave <- rowMeans(bulk_bin)
  return(sum(bulk_bin_ave > 0.65)/length(common.genes))
  
  
})
genes_expressed_sc_4_list <- sapply(colnames(strict), function(cell){
  rownames(strict[strict[,cell] > 0,])
})
strict_bulk_detection <- sapply(colnames(aggr_raw_GeTMM), function(cell){
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]
  common.genes <- intersect(rownames(bulk_raw_GeTMM), genes_expressed_sc_4_list[[cell]])
  
  bulk <- bulk_raw_GeTMM[common.genes, samples_cell_type %in% c(cell)]
  
  bulk_bin <- bulk > 5
  bulk_bin_ave <- rowMeans(bulk_bin)
  return(sum(bulk_bin_ave > 0.65)/length(common.genes))
  
  
})

data.frame(row.names = names(liberal_bulk_detection),
           cell = names(liberal_bulk_detection),
           #liberal = liberal_bulk_detection,
           medium = medium_bulk_detection,
           #conseravtive = conservative_bulk_detection,
           strict = strict_bulk_detection) %>% 
  reshape::melt(.) %>%
  ggplot() + geom_col(aes(x = cell, y = value, fill = variable), position = 'dodge')




genes_expressed_sc_2_list <- sapply(colnames(medium), function(cell){
  rownames(medium[medium[,cell] > 0,])
})




strict_TPM <- cengenDataSC::cengen_TPM_bulk[rownames(strict),colnames(strict)] * strict


strict_bulk_corr_pairwise <- pbsapply(colnames(aggr_raw_GeTMM), function(cell){

  
  #common.genes <- intersect(rownames(aggr_raw_GeTMM), rownames(strict_TPM))
  
  all_corr <- sapply(colnames(aggr_raw_GeTMM), function(cell2){
    
    common.genes <- intersect(rownames(aggr_raw_GeTMM), genes_expressed_sc_4_list[[cell2]])
    
    
    Bulk <- aggr_raw_GeTMM[common.genes,cell]
    SC <- strict_TPM[common.genes,cell2]
    
    return(cor(Bulk, SC, method = 'spearman'))
  })
  
  return(all_corr)
  
})


col_fun_ <- circlize::colorRamp2(colors = c('white', '#Cc0202'), breaks = c(min(unlist(strict_bulk_corr_pairwise)),
                                                                       max(unlist(strict_bulk_corr_pairwise))))

order(sc_size[])
sc_size['DD'] <- sc_size['VD_DD']
sc_size['VD'] <- sc_size['VD_DD']


col_dend = as.dendrogram(hclust(dist(t(strict_bulk_corr_pairwise))))


Heatmap(strict_bulk_corr_pairwise, cluster_rows = col_dend, cluster_columns = col_dend, row_names_side = 'left',
        col = col_fun_,
        show_row_dend = F, show_column_dend = F, name = ' ',
        row_title = 'Single Cell', column_title = 'Bulk', column_title_side = 'bottom',
        column_title_gp = gpar(fontsize = 30, fontface = "bold"), 
        row_title_gp = gpar(fontsize = 30, fontface = "bold"),
        column_names_gp = gpar(fontsize = 15, fontface = "bold"), 
        row_names_gp = gpar(fontsize = 15, fontface = "bold"))


