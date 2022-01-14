
### script to generate Figure 1 panels D, E, and F

### load the datasets, subset to ground truth genes and common neuron classes, and then calculate the TPR, FPR, and FDR for neuronal genes
### also calculate the FPR for non-neuronal genes for each dataset, plot the curves, and save.

library(tidyverse)
library(dplyr)
library(reshape)

library(wbData)
library(EnvStats)
library(pbapply)
library(bayestestR)

library(ggplot2)
library(ggrastr)
library(patchwork)


# ~~ functions ----

get_tpr <- function(expression, truth, threshold, na.rm = TRUE){
  # True Positive Rate, aka sensitivity, aka recall
  # TPR = TP/(TP+FN) = TP/P
  bin <- expression >= threshold
  return(sum(bin * truth)/sum(truth))
}
get_fpr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  bin <- expression >= threshold
  return(sum(bin * (!truth))/sum(!(truth)))
}
get_fdr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Discovery Rate
  # FDR = FP/(FP+TP) = 1 - PPV
  bin <- expression >= threshold
  fdr <- sum(bin * (!truth))/(sum(bin*(!truth)) + sum(bin*truth))
  if(is.nan(fdr))
    fdr <- 0
  return(fdr)
}

### ~~ function for integrating ----




# ~~ Ground truth ----

all_genes_bulk <- read.csv('bulk_all_ground_truth_042021.csv', row.names = 1)

set.seed(12345)


#~~ Single cell ----
CeNGEN_TPM <- read.table('CeNGEN_TPM_080421.tsv')
CeNGEN_TPM <- CeNGEN_TPM[order(rownames(CeNGEN_TPM)), order(colnames(CeNGEN_TPM))]



bulk_data <- read.table('bsn9_bulk_counts_113021.tsv', sep = '\t')

## if present, remove outliers and cell types with only 1 replicate
bulk_load$RICr133 <- NULL
bulk_load$PVMr122 <- NULL
bulk_load$ADFr99 <- NULL
bulk_load$M4r117 <- NULL



## if present, remove RNAs that are not being considered for technical reasons.
genes_keep <- ws281 %>% dplyr::filter(biotype != 'rRNA_gene' & biotype != 'miRNA_gene' & 
                                        biotype != 'piRNA_gene' & biotype != 'piRNA_gene' &
                                        biotype != 'transposable_element_gene' & biotype != 'antisense_lncRNA_gene' &
                                        biotype != 'scRNA_gene' & biotype != 'gene') %>% rownames()



bulk_data <- bulk_data[genes_keep,]


## load gene metadata table
bulk_meta <- read.table('bsn9_CeNGEN_bulk_geneLevel_metadata.tsv', sep = '\t')

#### normalize to gene length
bulk_data_pk <- (bulk_data/bulk_meta[rownames(bulk_data),'Length'])*1000


### load average integrated dataset
average_integration_GeTMM <- read.table('Average_integrated_TMM_counts_lengthNormalized_111521.tsv', sep = '\t')


### TMM normalize the unadjusted bulk counts
bulk_raw_GeTMM <- DGEList(counts = bulk_data_pk)
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM, method = 'TMM')
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)





### only consider cell types with samples in the integrated dataset
bulk_samples <- colnames(average_integration_pk)
bulk_sample_cell_types <- str_split_fixed(colnames(average_integration_pk), 'r', 2)[,1]
cell_types <- unique(bulk_sample_cell_types)
cell_types_wReps <- cell_types[base::table(bulk_sample_cell_types) > 1]
samples_wReps <- bulk_samples[bulk_sample_cell_types %in% cell_types_wReps]
samples_wReps_cell_types <- str_split_fixed(samples_wReps, 'r', 2)[,1]



# Make conformable gene list & neuron list
common.genes <- intersect(rownames(average_integration_GeTMM), rownames(CeNGEN_TPM))

bulk_neurs <- base::intersect(colnames(all_genes_bulk), cell_types_wReps)

###~~  average within cell types ----

aggr_raw_GeTMM <- bulk_raw_GeTMM
colnames(aggr_raw_GeTMM) <-str_split_fixed(colnames(aggr_raw_GeTMM),"r",2)[,1]
aggr_raw_GeTMM <- data.frame(vapply(unique(colnames(aggr_raw_GeTMM)), function(x) 
  rowMeans(aggr_raw_GeTMM[,colnames(aggr_raw_GeTMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_GeTMM)) ))
dim(aggr_raw_GeTMM)

aggr_ave_integrant_Getmm <- average_integration_GeTMM[,samples_wReps]
colnames(aggr_ave_integrant_Getmm) <-str_split_fixed(colnames(aggr_ave_integrant_Getmm),"r",2)[,1]
aggr_ave_integrant_Getmm <- data.frame(vapply(unique(colnames(aggr_ave_integrant_Getmm)), function(x) 
  rowMeans(aggr_ave_integrant_Getmm[,colnames(aggr_ave_integrant_Getmm)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_ave_integrant_Getmm)) ))








## ~~ get gene and neuron lists to use ----

neurons <- intersect(colnames(aggr_ave_integrant_Getmm), colnames(all_genes_bulk))
genes <- intersect(rownames(aggr_ave_integrant_Getmm), rownames(all_genes_bulk))


## ~~ subset data to just genes and neurons to use
aggr_ave_integrant_Getmm_plot <- aggr_ave_integrant_Getmm[genes, neurons]
aggr_raw_GeTMM_plot <- aggr_raw_GeTMM[genes, neurons]
CeNGEN_TPM_plot <- CeNGEN_TPM[genes, neurons]
all_gt <- all_genes_bulk[genes, neurons]



###~~ calc TPR FPR FDR ----



diags_CeNGEN_TPM_plot <- tibble(threshold = c(0,2**seq(-4,12,0.1)),
                                TPR = map_dbl(threshold, ~get_tpr(CeNGEN_TPM_plot, all_gt, .x)),
                                FPR = map_dbl(threshold, ~get_fpr(CeNGEN_TPM_plot, all_gt, .x)),
                                FDR = map_dbl(threshold, ~get_fdr(CeNGEN_TPM_plot, all_gt, .x)),
                                counts = "CeNGEN_TPM_plot")


diags_aggr_ave_integrant_Getmm_plot <- tibble(threshold = c(0,2**seq(-4,12,0.1)),
                                              TPR = map_dbl(threshold, ~get_tpr(aggr_ave_integrant_Getmm_plot, all_gt, .x)),
                                              FPR = map_dbl(threshold, ~get_fpr(aggr_ave_integrant_Getmm_plot, all_gt, .x)),
                                              FDR = map_dbl(threshold, ~get_fdr(aggr_ave_integrant_Getmm_plot, all_gt, .x)),
                                              counts = "aggr_ave_integrant_Getmm_plot")


diags_aggr_raw_GeTMM_plot <- tibble(threshold = c(0,2**seq(-4,12,0.1)),
                                    TPR = map_dbl(threshold, ~get_tpr(aggr_raw_GeTMM_plot, all_gt, .x)),
                                    FPR = map_dbl(threshold, ~get_fpr(aggr_raw_GeTMM_plot, all_gt, .x)),
                                    FDR = map_dbl(threshold, ~get_fdr(aggr_raw_GeTMM_plot, all_gt, .x)),
                                    
                                    counts = "aggr_raw_GeTMM_plot")



bind_rows(diags_CeNGEN_TPM_plot,
          diags_aggr_ave_integrant_Getmm_plot,
          diags_aggr_raw_GeTMM_plot) %>%
  ggplot(aes(x = FPR, y=TPR, color= counts)) +
  geom_point(size = 3) +
  theme_classic(base_size = 20) +
  theme(legend.position = '') +
  xlim(0,0.99)


bind_rows(diags_CeNGEN_TPM_plot[2:nrow(diags_CeNGEN_TPM_plot),],
          diags_aggr_ave_integrant_Getmm_plot[2:nrow(diags_aggr_ave_integrant_Getmm_plot),],
          diags_aggr_raw_GeTMM_plot[2:nrow(diags_aggr_raw_GeTMM_plot),]) %>%
  ggplot(aes(x = 1-FDR, y=TPR, color= counts)) +
  geom_point(size = 3) +
  theme_classic(base_size = 20) +
  theme(legend.position = '')


### calculate AUROC and AUPR

auc(diags_CeNGEN_TPM_plot$FPR, diags_CeNGEN_TPM_plot$TPR)
auc(diags_aggr_ave_integrant_Getmm_plot$FPR, diags_aggr_ave_integrant_Getmm_plot$TPR)
auc(diags_aggr_raw_GeTMM_plot$FPR, diags_aggr_raw_GeTMM_plot$TPR)

auc(1-diags_CeNGEN_TPM_plot$FDR, diags_CeNGEN_TPM_plot$TPR)
auc(1-diags_aggr_ave_integrant_Getmm_plot$FDR, diags_aggr_ave_integrant_Getmm_plot$TPR)
auc(1-diags_aggr_raw_GeTMM_plot$FDR, diags_aggr_raw_GeTMM_plot$TPR)







### load non-neuronal ground truth genes

likely_non_neuronal <- read.table('simpleMine_non_neuronal_genes_101521.tsv', sep = '\t',header = T)
likely_non_neuronal <- likely_non_neuronal[likely_non_neuronal$Putative.match.=='Yes',]

likely_non_neuronal_gt <- data.frame(row.names = likely_non_neuronal$WormBase.Gene.ID, 
                                     matrix(0, ncol=ncol(all_genes_bulk), nrow=nrow(likely_non_neuronal)))
colnames(likely_non_neuronal_gt) <- colnames(all_genes_bulk)

nn_genes <- intersect(rownames(likely_non_neuronal_gt), common.genes)
nn_genes <- intersect(nn_genes, rownames(prop_by_type_adjusted))

likely_non_neuronal_gt_cut <- likely_non_neuronal_gt[nn_genes, neurons]




### subset datasets to non-neuronal genes

CeNGEN_TPM_nnplot <- CeNGEN_TPM[nn_genes, neurons]
aggr_ave_integrant_Getmm_nnplot <- aggr_ave_integrant_Getmm[nn_genes, neurons]
aggr_raw_GeTMM_nnplot <- aggr_raw_GeTMM[nn_genes, neurons]



diags_CeNGEN_TPM_nnplot <- tibble(threshold = c(0,2**seq(-3,10,0.1)),
                                  FPR = map_dbl(threshold, ~get_fpr(CeNGEN_TPM_nnplot, likely_non_neuronal_gt_cut, .x)),
                                  counts = "CeNGEN_TPM_nnplot")

diags_aggr_ave_integrant_Getmm_nnplot <- tibble(threshold = c(0,2**seq(-3,10,0.1)),
                                                FPR = map_dbl(threshold, ~get_fpr(aggr_ave_integrant_Getmm_nnplot, likely_non_neuronal_gt_cut, .x)),
                                                counts = "aggr_ave_integrant_Getmm_nnplot")

diags_aggr_raw_GeTMM_nnplot <- tibble(threshold = c(0,2**seq(-3,10,0.1)),
                                      FPR = map_dbl(threshold, ~get_fpr(aggr_raw_GeTMM_nnplot, likely_non_neuronal_gt_cut, .x)),
                                      counts = "aggr_raw_GeTMM_nnplot")


#### plot the ROC, PR, and non-neuron FPR curves

bind_rows(diags_CeNGEN_TPM_plot,
          diags_aggr_ave_integrant_Getmm_plot,
          diags_aggr_raw_GeTMM_plot) %>%
  ggplot(aes(x = FPR, y=TPR, color= counts)) +
  geom_point_rast(size = 3) +
  theme_classic(base_size = 35) +
  theme(legend.position = '') +
  xlim(0,0.99) +
  
  bind_rows(diags_CeNGEN_TPM_plot[2:nrow(diags_CeNGEN_TPM_plot),],
            diags_aggr_ave_integrant_Getmm_plot[2:nrow(diags_aggr_ave_integrant_Getmm_plot),],
            diags_aggr_raw_GeTMM_plot[2:nrow(diags_aggr_raw_GeTMM_plot),]) %>%
  ggplot(aes(x = 1-FDR, y=TPR, color= counts)) +
  geom_point_rast(size = 3) +
  xlab('Precision') +
  theme_classic(base_size = 35) +
  theme(legend.position = '') +
  
  
  bind_rows(diags_CeNGEN_TPM_nnplot,
            diags_aggr_raw_GeTMM_nnplot,
            diags_aggr_ave_integrant_Getmm_nnplot) %>%
  ggplot(aes(x = threshold, y=FPR, color= counts)) +
  geom_point_rast(size = 3) + xlim(0, 100) +
  theme_classic(base_size = 35) +
  theme(legend.position = '')

ggsave('ROC_PR_nnFPR_curves_113021.pdf', width = 30, height = 14)


