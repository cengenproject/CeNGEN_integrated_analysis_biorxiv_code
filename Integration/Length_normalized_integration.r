## if needed, install libraries
install.packages('stringr')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("AlexWeinreb/wbData")

### load libraries
library(stringr)
library(wbData)
library(edgeR)

## function for integration

### inputs:
# bulk genes x samples matrix. Colnames should be: cell_type#sampleID, where # is in place of some string to split on. default currently is lower case r
# single cell pseudobulk genes x samples matrix Colnames should be: cell_type#sampleID, where # is in place of some string to split on. default currently is two underscores
# pseudocount, no default given, I use 0.1 to reduce distortion of non-zero counts

integrate_geometricMean_biorep <- function(bulk_replicates, single_cell_replicates, pseudocount, bulk_sep = 'r', single_cell_sep = '__'){

  common.genes <- intersect(rownames(bulk_replicates), rownames(single_cell_replicates))

  bulk_replicates <- bulk_replicates[common.genes,]
  single_cell_replicates <- single_cell_replicates[common.genes,]

  bulk_cell_types <- str_split_fixed(colnames(bulk_replicates), bulk_sep, 2)[,1]
  cells_bio_rep <- str_split_fixed(colnames(single_cell_replicates), single_cell_sep, 2)[,1]

  index_list_sc <- numeric()
  index_list_bulk <- numeric()
  for(g in unique(bulk_cell_types)){
    rep_totals <- c(sum(bulk_cell_types==g), sum(cells_bio_rep==g))### count how many replicates are in each cell type for each data source
    rep_min <- min(rep_totals) ## select the smaller number, this will be the number of samples in the final integrated matrixes
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

  single_cell_replicates_match <- single_cell_replicates[,index_list_sc] ### normalize the Single Cell data to match the bulk data
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


## load single cell aggregates, precomputed. Each is the sum of counts for cells in a cluster-replicate (ex: all AFD cells that came from replicate cho.1_2)
## limited to cluster-replicates with >= 10 cells

CeNGEN_aggr_bio_reps <- read.table('~/Bioinformatics/single_cell_data/CeNGEN_aggr_bio_reps_101821.tsv')

print('loaded single cell reps')



CeNGEN_aggr_bio_reps <- CeNGEN_aggr_bio_reps[,order(colnames(CeNGEN_aggr_bio_reps))]


### expand to VD and DD, with the caveat that we then can't do differential expression comparing VD and DD using this dataset
cell_labels <- str_split_fixed(colnames(CeNGEN_aggr_bio_reps), '__',2)[,1]
CeNGEN_aggr_bio_reps_DD.temp <- CeNGEN_aggr_bio_reps[,cell_labels == 'VD_DD']
CeNGEN_aggr_bio_reps_VD.temp <- CeNGEN_aggr_bio_reps[,cell_labels == 'VD_DD']
colnames(CeNGEN_aggr_bio_reps_DD.temp) <- paste0('DD__', str_split_fixed(colnames(CeNGEN_aggr_bio_reps_DD.temp), '__', 2)[,2])
colnames(CeNGEN_aggr_bio_reps_VD.temp) <- paste0('VD__', str_split_fixed(colnames(CeNGEN_aggr_bio_reps_VD.temp), '__', 2)[,2])


CeNGEN_aggr_bio_reps.test <- cbind(CeNGEN_aggr_bio_reps, CeNGEN_aggr_bio_reps_DD.temp, CeNGEN_aggr_bio_reps_VD.temp)
CeNGEN_aggr_bio_reps <- CeNGEN_aggr_bio_reps.test[,order(colnames(CeNGEN_aggr_bio_reps.test))]
cell_labels <- str_split_fixed(colnames(CeNGEN_aggr_bio_reps), '__',2)[,1]


bulk_cell_types <- str_split_fixed(colnames(bulk_data), 'r', 2)[,1]


### get the bulk metadata, to normalize genes by length prior to integration, if desired.
### not recommended if performing differential expression! This is primarily for increased accuracy in calling expression.

bulk_meta <- read.table('~/Bioinformatics/bsn5/bsn5_bulk_counts_metadata.tsv')
bulk_data_pk <- (bulk_data/bulk_meta[rownames(bulk_data), 'Length']) * 1000


## which cell types have > 1 replicate? Performed for both bulk and pseudobulk replicates

samples <- colnames(bulk_data)
samples_cell_types <- do.call(rbind, strsplit(colnames(bulk_data), 'r', 2))[,1]
cell_types <- unique(samples_cell_types)
cell_types_wReps <- cell_types[base::table(samples_cell_types) > 1]
samples_wReps <- samples[samples_cell_types %in% cell_types_wReps]
samples_wReps_cell_types <- do.call(rbind, strsplit(samples_wReps, 'r', 2))[,1]

samples_SC <- colnames(CeNGEN_aggr_bio_reps)
samples_cell_types_SC <- do.call(rbind, strsplit(colnames(CeNGEN_aggr_bio_reps), '__', 2))[,1]
cell_types_SC <- unique(samples_cell_types_SC)
cell_types_wReps_SC <- cell_types_SC[base::table(samples_cell_types_SC) > 1]
samples_wReps_SC <- samples_SC[samples_cell_types_SC %in% cell_types_wReps_SC]
samples_wReps_cell_types_SC <- do.call(rbind, strsplit(samples_wReps_SC, '__', 2))[,1]



### integrate them, and if desired, save them

many_integrations_pk <- pblapply(seq(101,150,1), function(v){
  set.seed(v)
  test_int <- integrate_geometricMean_biorep(bulk_data_pk[,samples_wReps], CeNGEN_aggr_bio_reps[,samples_wReps_SC],
                          pseudocount = 0.1, bulk_sep = 'r', single_cell_sep = '__')
  return(test_int) })
print('integrated samples')

colnames(many_integrations_pk[[1]])



#### TMM normalize the integrated samples
many_GeTMMs <- pblapply(many_integrations_pk, function(integrant){
  GeTMM <- DGEList(integrant)
  GeTMM <- calcNormFactors(GeTMM)
  GeTMM <- cpm(GeTMM, normalized.lib.sizes = T)
  GeTMM[GeTMM < 0.0000000000001] <- 0
  return(GeTMM)
})
many_GeTMMs[[1]][1:10,1:10]
bulk_raw_GeTMM[1:10,1:10]


### get the average normalized samples
ave_integrated_GeTMM <- Reduce('+', many_GeTMMs)/length(many_GeTMMs)
ave_integrated_GeTMM <- ave_integrated_GeTMM[order(rownames(ave_integrated_GeTMM)),
                                                             order(colnames(ave_integrated_GeTMM))]
                                                             
### get the average cell_type profile                                                             
aggr_ave_integrated_GeTMM <- ave_integrated_GeTMM
colnames(aggr_ave_integrated_GeTMM) <-str_split_fixed(colnames(aggr_ave_integrated_GeTMM),"r",2)[,1]
aggr_ave_integrated_GeTMM <- data.frame(vapply(unique(colnames(aggr_ave_integrated_GeTMM)), function(x) 
  rowMeans(aggr_ave_integrated_GeTMM[,colnames(aggr_ave_integrated_GeTMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_ave_integrated_GeTMM)) ))


### save the integrated matrixes if desired. Here written as separate csvs. Easy to change if you want to save as 1 RDS file.
for(g in seq(1,50,1)){
  write.csv(many_integrations1[[g]], paste0('integrated_matrix_', g, '.csv'))
}
