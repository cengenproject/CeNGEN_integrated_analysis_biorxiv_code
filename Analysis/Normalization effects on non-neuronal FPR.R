###

### libraries ----
library(ggplot2)
library(ggrastr)
library(edgeR)
library(bayestestR)

## function ----

get_fpr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  bin <- expression >= threshold
  return(sum(bin * (!truth))/sum(!(truth)))
}

#### load ground truth matrix ----
nonNeuronal_GT <- read.csv('bulk_nonNeuronal_ground_truth_040122.csv', row.names = 1)

#### load bulk data ----
bulk_data <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.tsv')

bulk_meta <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_geneLevel_metadata.tsv')

bulk_meta <- bulk_meta[rownames(bulk_data),]

bulk_data_pk <- (bulk_data/bulk_meta$Length) * 1000

#### comparing TMM and GeTMM normalized counts for non-neuronal FPR ----

bulk_raw_GeTMM <- DGEList(counts = bulk_data_pk)
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM, method = 'TMM')
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)

bulk_raw_TMM <- DGEList(counts = bulk_data)
bulk_raw_TMM <- calcNormFactors(bulk_raw_TMM, method = 'TMM')
bulk_raw_TMM <- cpm(bulk_raw_TMM, normalized.lib.sizes = T)

### get the average profiles per cell type ----
aggr_raw_GeTMM <- bulk_raw_GeTMM
colnames(aggr_raw_GeTMM) <-str_split_fixed(colnames(aggr_raw_GeTMM),"r",2)[,1]
aggr_raw_GeTMM <- data.frame(vapply(unique(colnames(aggr_raw_GeTMM)), function(x) 
  rowMeans(aggr_raw_GeTMM[,colnames(aggr_raw_GeTMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_GeTMM)) ))
dim(aggr_raw_GeTMM)


aggr_raw_TMM <- bulk_raw_TMM
colnames(aggr_raw_TMM) <-str_split_fixed(colnames(aggr_raw_TMM),"r",2)[,1]
aggr_raw_TMM <- data.frame(vapply(unique(colnames(aggr_raw_TMM)), function(x) 
  rowMeans(aggr_raw_TMM[,colnames(aggr_raw_TMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_TMM)) ))




### cut down to just neurons with samples
colnames(nonNeuronal_GT) <- colnames(aggr_raw_GeTMM)

nn_genes <- intersect(rownames(nonNeuronal_GT), aggr_raw_GeTMM)


nonNeuronal_GT <- nonNeuronal_GT[nn_genes,]
aggr_raw_GeTMM_nnplot <- aggr_raw_GeTMM[nn_genes,]
aggr_raw_TMM_nnplot <- aggr_raw_TMM[nn_genes,]


### get FPR for both datasets across a wide range of thresholds

diags_aggr_raw_TMM_nnplot <- tibble(threshold = c(0,2**seq(-3,10,0.1)),
                                                FPR = map_dbl(threshold, ~get_fpr(aggr_raw_TMM_nnplot, nonNeuronal_GT, .x)),
                                                counts = "aggr_raw_TMM_nnplot")

diags_aggr_raw_GeTMM_nnplot <- tibble(threshold = c(0,2**seq(-3,10,0.1)),
                                      FPR = map_dbl(threshold, ~get_fpr(aggr_raw_GeTMM_nnplot, nonNeuronal_GT, .x)),
                                      counts = "aggr_raw_GeTMM_nnplot")

## plot
bind_rows(diags_aggr_raw_TMM_nnplot,
          diags_aggr_raw_GeTMM_nnplot) %>%
  ggplot(aes(x = log10(threshold+1), y=FPR, color= counts)) +
  geom_point_rast(size = 3) +
  scale_x_continuous(breaks = log10(c(1,11,101,1001, 10001)),
                     labels = c('0', '10', '100', '1,000', '10k'))+
  
  theme_classic(base_size = 35) +
  theme(legend.position = '') +
  xlab('Threshold log10(x+1)')

ggsave('Normalized_vs_unNormalized_bulk_nnFPR_011928.pdf')


## calculate FPR AUC
auc(log10(diags_aggr_raw_TMM_nnplot$threshold+1)/log10(max(diags_aggr_raw_TMM_nnplot$threshold+1)), diags_aggr_raw_TMM_nnplot$FPR)
auc(log10(diags_aggr_raw_GeTMM_nnplot$threshold+1)/log10(max(diags_aggr_raw_GeTMM_nnplot$threshold+1)), diags_aggr_raw_GeTMM_nnplot$FPR)



