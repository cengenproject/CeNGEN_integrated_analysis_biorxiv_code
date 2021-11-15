### expression of protein coding genes not detected in any single cell clusters


### datasets needed, and assumed already loaded: ----

# CeNGEN_TPM: bulk annotation aggregates for single cell data
# Aggr_raw_GeTMM: aggregated raw bulk GeTMM values, arithmetic mean per cell type
# GeTMM_cor_s_contaminant: genes x contaminants spearman correlations
# ws277: wormbase gene metadata, from wbData package
# nn_genes: a list of non-neuronal genes, curated from wormbase using simpleMine


### load Alexis single cell thresholds

ws277_prc <- ws277[ws277$biotype=='protein_coding_gene',]

sc_gene_thresholds <- readRDS('~/211028_genes_categorized_by_pattern.rds')

nondetected_threshold_genes <- sc_gene_thresholds[['nondetected']]
nondetected_threshold_genes <- intersect(nondetected_threshold_genes, rownames(ws277_prc))


not_detected_Seurat <- setdiff(rownames(aggr_raw_GeTMM), rownames(CeNGEN_TPM))
not_detected_Seurat <- intersect(not_detected_Seurat, rownames(ws277_prc))
length(not_detected_Seurat)


### subset to the genes of interest given either criteria and combine the matrices ----
aggr_not_detected_Seurat <- aggr_raw_GeTMM[not_detected_Seurat,]
aggr_nondetected_threshold_genes <- aggr_raw_GeTMM[nondetected_threshold_genes,]

aggr_nondetected_all <- rbind(aggr_nondetected_threshold_genes, aggr_not_detected_Seurat)
aggr_nondetected_all <- na.omit(aggr_nondetected_all)
dim(aggr_nondetected_all)


### subset the correlation to contaminants ----
GeTMM_cor_s_contaminant_prc <- GeTMM_cor_s_contaminant[intersect(rownames(GeTMM_cor_s_contaminant), rownames(ws277_prc)),]

GeTMM_cor_s_contaminant_prc_nondetect_all <- GeTMM_cor_s_contaminant_prc[intersect(rownames(GeTMM_cor_s_contaminant_prc), 
                                                                                   rownames(aggr_nondetected_all)),]


### plot the maximum correlation per gene
data.frame( max_corr = apply(GeTMM_cor_s_contaminant_prc_nondetect_all, 1, max)) |> ggplot() + geom_density(aes(x = max_corr))


### keep just the genes with a low correlation
low_corr_genes_keep <- rownames(GeTMM_cor_s_contaminant_prc_nondetect_all)[apply(GeTMM_cor_s_contaminant_prc_nondetect_all, 1, max) < 0.3]



### Check that most non-neuronal genes are removed by thresholding on the correlation to contaminants

length(intersect(nn_genes, rownames(aggr_nondetected_all)))
low_expression_nn_genes2 <- intersect(nn_genes, low_corr_genes_keep)
length(low_expression_nn_genes2)


### set up ground truth
low_expression_likely_non_neuronal_gt_cut2 <- likely_non_neuronal_gt[low_expression_nn_genes2, neurons]
aggr_raw_GeTMM_low2_nnplot <- aggr_raw_GeTMM[low_expression_nn_genes2, neurons]


### calculate FPR rate for the remaining non-neuronal genes
diags_aggr_raw_GeTMM_low2_nnplot <- tibble(threshold = c(0,2**seq(-4,12,0.1)),
                                           FPR = map_dbl(threshold, ~get_fpr(aggr_raw_GeTMM_low2_nnplot, 
                                                                             low_expression_likely_non_neuronal_gt_cut2, .x)),
                                           counts = "aggr_raw_GeTMM_low2_nnplot")

ggplot(diags_aggr_raw_GeTMM_low2_nnplot[diags_aggr_raw_GeTMM_low2_nnplot$threshold < 100,], aes(x = threshold, y = FPR)) + 
  geom_point() +xlim(0,100) + geom_hline(yintercept = 0) + theme_classic(base_size = 20) +
  xlab('GeTMM threshold') + ylab('Non-Neuronal FPR')



## find a threshold that excludes all known non-neuronal genes
diags_aggr_raw_GeTMM_low2_nnplot[diags_aggr_raw_GeTMM_low2_nnplot$FPR == 0,]


### calculate the new genes per cell type ----
new_genes <- sapply(colnames(aggr_nondetected_all), function(cell){
  sum(aggr_nondetected_all[low_corr_genes_keep,cell] > 16)
})
new_genes


## bar plot ----
data.frame(new_genes,
           cell_type = names(new_genes)) |>
  ggplot() +
  geom_col(aes(x = cell_type, y = new_genes)) +
  theme_classic(base_size = 20)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold')) +
  xlab('Cell Type') + ylab('"New" Genes Detected in Bulk\n(Never Detected in Single Cell)')
