### do pairwise differential expression on all bulk uncorrected samples

## load libraries ----
library(tidyverse)
library(dplyr)
library(wbData)
library(edgeR)
library(pbapply)
library(reshape)
library(stringr)
library(wbData)

print('loaded libraries')

bulk_data <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.tsv')

print('loaded bulk data')

#############
## load in reference information on gene classes
ref_genes <- data.frame(wbData::wb_load_gene_ids('281'))
rownames(ref_genes) <- ref_genes$gene_id



ref_genes_protein_coding <- ref_genes[ref_genes$biotype == 'protein_coding_gene',]



bulk_data <- bulk_data[rownames(bulk_data) %in% rownames(ref_genes_protein_coding),order(colnames(bulk_data))]

print('subsetted to protein coding genes')

cell_labels <- str_split_fixed(colnames(bulk_data), 'r',2)[,1]


## which cell types have > 1 replicate?

samples <- colnames(bulk_data)
samples_cell_types <- cell_labels
cell_types <- unique(samples_cell_types)
cell_types_wReps <- cell_types[base::table(samples_cell_types) > 1]
samples_wReps <- samples[samples_cell_types %in% cell_types_wReps]
samples_wReps_cell_types <- do.call(rbind, strsplit(samples_wReps, 'r', 2))[,1]



### running edgeR

bulk_edgeR <- edgeR::DGEList(counts = bulk_data[,samples_wReps], group = samples_wReps_cell_types)
bulk_edgeR <- calcNormFactors(bulk_edgeR, method = 'TMM')

print('calculated TMM')

group <- factor(bulk_edgeR$samples$group)
design_no_int <- model.matrix(~0+group, data=bulk_edgeR$samples)
colnames(design_no_int) <- levels(bulk_edgeR$samples$group)


bulk_edgeR <- estimateDisp(bulk_edgeR, design_no_int)
print('estimated dispersion')
bulk_fit <- glmQLFit(bulk_edgeR, design_no_int)
print('fit finished')


saveRDS(bulk_edgeR, 'bulk_edgeR_112821.rds')
saveRDS(bulk_fit, 'bulk_fit_112821.rds')

print('fit model')


### making contrasts

contrasts_names <- factor()
cells_shrink <- unique(bulk_fit$samples$group)

#design matrix
design_no_int <- model.matrix(~0+factor(bulk_fit$samples$group), data=bulk_fit$samples)
colnames(design_no_int) <- levels(bulk_fit$samples$group)

# all neuron-neuron pairs, just one per pair, not both directions
for(f in cells_shrink){
  cells_shrink <- cells_shrink[!(cells_shrink %in% f)]
  for(g in cells_shrink){
    contrasts_names <- c(contrasts_names, paste0(f,'-(',g,')'))
  }
}
names(contrasts_names) <- contrasts_names

contrasts_edgeR <- lapply(contrasts_names, function(v){
  cont <- limma::makeContrasts(contrasts = contrasts_names[[v]], levels = design_no_int)
})
names(contrasts_edgeR) <- contrasts_names
print('made contrasts')

### run qlf test
qlf_tables <- lapply(names(contrasts_edgeR), function(v){
  print(v)
  qlf <- glmQLFTest(glmfit = bulk_fit, contrast = contrasts_edgeR[[v]])
  return(qlf$table)
})
names(qlf_tables) <- names(contrasts_edgeR)
print('qlf done')
saveRDS(qlf_tables, 'bulk_data_qlf_tables_bsn9_112821.rds')
