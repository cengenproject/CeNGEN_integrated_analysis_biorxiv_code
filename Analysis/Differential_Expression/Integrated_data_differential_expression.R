
library(tidyverse)
library(dplyr)
library(wbData)
library(edgeR)
library(pbapply)
library(reshape)
library(stringr)
library(wbData)


#### integration function ----

integrate_geometricMean_biorep <- function(bulk_replicates, single_cell_replicates, pseudocount, bulk_sep = 'r', single_cell_sep = '__'){

  common.genes <- intersect(rownames(bulk_replicates), rownames(single_cell_replicates))

  bulk_replicates <- bulk_replicates[common.genes,]
  single_cell_replicates <- single_cell_replicates[common.genes,]

  bulk_cell_types <- str_split_fixed(colnames(bulk_replicates), bulk_sep, 2)[,1]
  cells_bio_rep <- str_split_fixed(colnames(single_cell_replicates), single_cell_sep, 2)[,1]

  index_list_sc <- numeric()
  index_list_bulk <- numeric()
  for(g in unique(bulk_cell_types)){
    rep_totals <- c(sum(bulk_cell_types==g), sum(cells_bio_rep==g))
    rep_min <- min(rep_totals)
    #print(g)
    #print(sum(bulk_cell_types==g))
    #print(rep_min)
    indices_sc <- sample(which(cells_bio_rep==g), size = rep_min)
    index_list_sc <- c(index_list_sc, indices_sc)
    indices_bulk <- sample(which(bulk_cell_types==g), size = rep_min)
    index_list_bulk <- c(index_list_bulk, indices_bulk)

  }
  #print(index_list_sc)
  #print(index_list_bulk)

  #print('set indexes')
  bulk_replicates_match <- bulk_replicates[,index_list_bulk]

  single_cell_replicates_match <- single_cell_replicates[,index_list_sc]
  single_cell_replicates_match <- sweep(single_cell_replicates_match, 2, colSums(single_cell_replicates_match, na.rm = T), '/') * colSums(bulk_replicates_match, na.rm = T)

  #bulk_replicates_match <- bulk_replicates_match[,order(colnames(bulk_replicates_match))]
  #single_cell_replicates_match <- single_cell_replicates_match[,order(colnames(single_cell_replicates_match))]
  #print(colnames(bulk_replicates_match))
  #print(colnames(single_cell_replicates_match))

  sample_level_integration <- exp( ( log(bulk_replicates_match + pseudocount) +
                                       log(single_cell_replicates_match + pseudocount) ) /2 ) - pseudocount

  sample_level_integration <- sample_level_integration[order(rownames(sample_level_integration)),
                                                       order(colnames(sample_level_integration))]

  return(sample_level_integration)
}




############################ ----------------- #################################

CeNGEN_aggr_bio_reps <- read.table('CeNGEN_aggr_bio_reps_101821.tsv')

cell_labels <- str_split_fixed(colnames(CeNGEN_aggr_bio_reps), '__',2)[,1]

#############



### remove rRNA, miRNA, and piRNA genes

## load in reference information on gene classes
ref_genes <- data.frame(wbData::wb_load_gene_ids('281'))
rownames(ref_genes) <- ref_genes$gene_id

### load in bulk data ----
bulk_data <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.tsv')


ref_genes_protein_coding <- ref_genes[ref_genes$biotype == 'protein_coding_gene',]



CeNGEN_aggr_bio_reps <- CeNGEN_aggr_bio_reps[,order(colnames(CeNGEN_aggr_bio_reps))]

cell_labels <- str_split_fixed(colnames(CeNGEN_aggr_bio_reps), '__',2)[,1]
CeNGEN_aggr_bio_reps_DD.temp <- CeNGEN_aggr_bio_reps[,cell_labels == 'VD_DD']
CeNGEN_aggr_bio_reps_VD.temp <- CeNGEN_aggr_bio_reps[,cell_labels == 'VD_DD']
colnames(CeNGEN_aggr_bio_reps_DD.temp) <- paste0('DD__', str_split_fixed(colnames(CeNGEN_aggr_bio_reps_DD.temp), '__', 2)[,2])
colnames(CeNGEN_aggr_bio_reps_VD.temp) <- paste0('VD__', str_split_fixed(colnames(CeNGEN_aggr_bio_reps_VD.temp), '__', 2)[,2])

CeNGEN_aggr_bio_reps.test <- cbind(CeNGEN_aggr_bio_reps, CeNGEN_aggr_bio_reps_DD.temp, CeNGEN_aggr_bio_reps_VD.temp)
CeNGEN_aggr_bio_reps <- CeNGEN_aggr_bio_reps[,order(colnames(CeNGEN_aggr_bio_reps))]
cell_labels <- str_split_fixed(colnames(CeNGEN_aggr_bio_reps), '__',2)[,1]

bulk_cell_types <- str_split_fixed(colnames(bulk_data), 'r', 2)[,1]


common.genes <- intersect(rownames(bulk_data), rownames(CeNGEN_aggr_bio_reps))
common.genes <- intersect(common.genes, ref_genes_protein_coding)


bulk_data <- bulk_data[common.genes,]
CeNGEN_aggr_bio_reps <- CeNGEN_aggr_bio_reps[common.genes,]


## which cell types have > 1 replicate?

samples <- colnames(bulk_data)
samples_cell_types <- do.call(rbind, strsplit(colnames(bulk_data), 'r', 2))[,1]
cell_types <- unique(samples_cell_types)
cell_types_wReps <- cell_types[base::table(samples_cell_types) > 1]
samples_wReps <- samples[samples_cell_types %in% cell_types_wReps]
samples_wReps_cell_types <- do.call(rbind, strsplit(samples_wReps, 'r', 2))[,1]

group <- factor(samples_wReps_cell_types)
group
design_no_int <- model.matrix(~0+group, data=data.frame(samples_wReps))
colnames(design_no_int) <- cell_types_wReps
rownames(design_no_int) <- samples_wReps



contrasts_names <- factor()
cells_shrink <- cell_types_wReps
for(f in cells_shrink){
  cells_shrink <- cells_shrink[!(cells_shrink %in% f)]
  for(g in cells_shrink){
    contrasts_names <- c(contrasts_names, paste0(f,'-(',g,')'))
  }
}
names(contrasts_names) <- contrasts_names

contrasts_edgeR <- lapply(contrasts_names, function(v){ cont <-
  makeContrasts(contrasts = contrasts_names[[v]], levels = design_no_int)
})

print('set up contrasts')


many_integrations <- pblapply(seq(101,151,1), function(v){

  set.seed(v)

  test_int <- integrate_geometricMean_biorep(bulk_data[,samples_wReps], CeNGEN_aggr_bio_reps, pseudocount = 0.1, bulk_sep = 'r', single_cell_sep = '__')
  return(test_int) })
print('integrated samples')



int_cell_types <- unique(str_split_fixed(colnames(many_integrations[[1]]), 'r', 2)[,1])
contrasts_names <- factor()

cells_shrink <- int_cell_types
for(f in cells_shrink){
  cells_shrink <- cells_shrink[!(cells_shrink %in% f)]
  for(g in cells_shrink){
    contrasts_names <- c(contrasts_names, paste0(f,'-(',g,')'))
  }
}
names(contrasts_names) <- contrasts_names

contrasts_edgeR <- lapply(contrasts_names, function(v){ cont <-
  makeContrasts(contrasts = contrasts_names[[v]], levels = design_no_int)
})

print('set up contrasts')


listing <- 1
many_edgeRs <- pblapply(many_integrations, function(integrated){

  samples <- colnames(integrated)
  samples_cell_types <- do.call(rbind, strsplit(colnames(integrated), 'r', 2))[,1]
  cell_types <- unique(samples_cell_types)
  cell_types_wReps <- cell_types[base::table(samples_cell_types) > 1]
  samples_wReps <- samples[samples_cell_types %in% cell_types_wReps]
  samples_wReps_cell_types <- do.call(rbind, strsplit(samples_wReps, 'r', 2))[,1]

  cell_labels <- str_split_fixed(colnames(integrated[,samples_wReps]), 'r', 2)[,1]
  dge <- DGEList(integrated[,samples_wReps], group = cell_labels)

  group <- factor(samples_wReps_cell_types)
  group
  design_no_int <- model.matrix(~0+group, data=data.frame(samples_wReps))
  colnames(design_no_int) <- cell_types_wReps
  rownames(design_no_int) <- samples_wReps

  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = design_no_int)

  print(listing)
  listing <<- listing + 1

  return(dge)  })

print('calculated dispersions')


cell_types_wReps <- colnames(many_edgeRs[[1]]$design)


contrasts_names <- factor()
cells_shrink <- cell_types_wReps
for(f in cells_shrink){
  cells_shrink <- cells_shrink[!(cells_shrink %in% f)]
  for(g in cells_shrink){
    contrasts_names <- c(contrasts_names, paste0(f,'-(',g,')'))
  }
}
names(contrasts_names) <- contrasts_names

contrasts_edgeR <- lapply(contrasts_names, function(v){ cont <-
  makeContrasts(contrasts = contrasts_names[[v]], levels = cell_types_wReps)
})




many_fits <- pblapply(many_edgeRs, function(dge){

  dge_fit <- glmQLFit(dge, dge$design)

  return(dge_fit)  })

print('calculated fits')



many_p_values_and_logFCs <- lapply(contrasts_edgeR, function(contrast){
  qlfs_temp <- pblapply(many_fits, function(fit){
    qlf <- glmQLFTest(glmfit=fit, contrast=contrast)
    return(qlf$table[,c('logFC', 'PValue')])
  })
  PValues <- lapply(qlfs_temp, function(tab){
    return(tab$PValue)
  })
  PValues <- do.call(cbind, PValues)
  logFCs <- lapply(qlfs_temp, function(tab){
    return(tab$logFC)
  })
  logFCs <- do.call(cbind, logFCs)

  return(list(PValues, logFCs))
})
names(many_p_values_and_logFCs) <- contrasts_edgeR

print('done')


saveRDS(many_p_values_and_logFCs, 'many_p_values_and_logFCs_5_112821.rds')
