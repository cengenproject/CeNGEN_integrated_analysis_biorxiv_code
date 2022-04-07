### detection of single cell genes in bulk
library(dendextend)
library(stringr)
library(ComplexHeatmap)
library(wbData)
library(edgeR)

### load single cell thresholded data with bulk level annotations (less resolution than single cell)
strict <- cengenDataSC::cengen_sc_4_bulk
strict[strict > 0] = 1


## get list of genes expressed in each cell type using single cell data
genes_expressed_sc_4_list <- sapply(colnames(strict), function(cell){
  rownames(strict[strict[,cell] > 0,])
})

## get thresholded TPM data
strict_TPM <- cengenDataSC::cengen_TPM_bulk[rownames(strict),colnames(strict)] * strict

## load bulk data and normalize to gene length
bulk_data <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.tsv', sep = '\t')

bulk_meta <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_geneLevel_metadata.tsv')

bulk_meta <- bulk_meta[rownames(bulk_data),]

bulk_data_pk <- (bulk_data/bulk_meta$Length) * 1000

## inter sample normalization using edgeR
bulk_raw_GeTMM <- DGEList(counts = bulk_data_pk)
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM, method = 'TMM')
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)


## get the average for each cell type
dim(bulk_raw_GeTMM)
aggr_raw_GeTMM <- bulk_raw_GeTMM
colnames(aggr_raw_GeTMM) <-str_split_fixed(colnames(aggr_raw_GeTMM),"r",2)[,1]
aggr_raw_GeTMM <- data.frame(vapply(unique(colnames(aggr_raw_GeTMM)), function(x) 
  rowMeans(aggr_raw_GeTMM[,colnames(aggr_raw_GeTMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_GeTMM)) ))


## calculate the Spearman correlation for each single cell type with each bulk cell type, using only the genes called expressed in the single cell data.
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

## plot heatmap
Heatmap(strict_bulk_corr_pairwise, cluster_rows = col_dend, cluster_columns = col_dend, row_names_side = 'left',
        col = col_fun_,
        show_row_dend = F, show_column_dend = F, name = ' ',
        row_title = 'Single Cell', column_title = 'Bulk', column_title_side = 'bottom',
        column_title_gp = gpar(fontsize = 30, fontface = "bold"), 
        row_title_gp = gpar(fontsize = 30, fontface = "bold"),
        column_names_gp = gpar(fontsize = 15, fontface = "bold"), 
        row_names_gp = gpar(fontsize = 15, fontface = "bold"))


