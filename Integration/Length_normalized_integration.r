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
library(pbapply)

## function for integration

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

  return(sample_level_integration)}


## load single cell aggregates

CeNGEN_aggr_bio_reps <- read.table('CeNGEN_aggr_bio_reps_101821.tsv')

CeNGEN_aggr_bio_reps <- CeNGEN_aggr_bio_reps[,order(colnames(CeNGEN_aggr_bio_reps))]


print('loaded single cell reps')

## load bulk data

bulk_data <- read.table('bsn9_CeNGEN_bulk_data.tsv', header = T)




## fix the single cell data to match bulk_data labels

cell_labels <- str_split_fixed(colnames(CeNGEN_aggr_bio_reps), '__',2)[,1]
CeNGEN_aggr_bio_reps_DD.temp <- CeNGEN_aggr_bio_reps[,cell_labels == 'VD_DD']
CeNGEN_aggr_bio_reps_VD.temp <- CeNGEN_aggr_bio_reps[,cell_labels == 'VD_DD']
colnames(CeNGEN_aggr_bio_reps_DD.temp) <- paste0('DD__', str_split_fixed(colnames(CeNGEN_aggr_bio_reps_DD.temp), '__', 2)[,2])
colnames(CeNGEN_aggr_bio_reps_VD.temp) <- paste0('VD__', str_split_fixed(colnames(CeNGEN_aggr_bio_reps_VD.temp), '__', 2)[,2])


CeNGEN_aggr_bio_reps.test <- cbind(CeNGEN_aggr_bio_reps, CeNGEN_aggr_bio_reps_DD.temp, CeNGEN_aggr_bio_reps_VD.temp)
CeNGEN_aggr_bio_reps <- CeNGEN_aggr_bio_reps.test[,order(colnames(CeNGEN_aggr_bio_reps.test))]
cell_labels <- str_split_fixed(colnames(CeNGEN_aggr_bio_reps), '__',2)[,1]



## get bulk data labels
bulk_cell_types <- str_split_fixed(colnames(bulk_data), 'r', 2)[,1]


## which cell types have > 1 replicate?

samples <- colnames(bulk_data)
samples_cell_types <- do.call(rbind, strsplit(colnames(bulk_data), 'r', 2))[,1]
cell_types <- unique(samples_cell_types)
cell_types_wReps <- cell_types[base::table(samples_cell_types) > 1]
samples_wReps <- samples[samples_cell_types %in% cell_types_wReps]
samples_wReps_cell_types <- do.call(rbind, strsplit(samples_wReps, 'r', 2))[,1]


sc_samples <- colnames(CeNGEN_aggr_bio_reps)
sc_samples_cell_types <- do.call(rbind, strsplit(colnames(CeNGEN_aggr_bio_reps), '__', 2))[,1]
sc_cell_types <- unique(sc_samples_cell_types)
sc_cell_types_wReps <- sc_cell_types[base::table(sc_samples_cell_types) > 1]
sc_samples_wReps <- sc_samples[sc_samples_cell_types %in% sc_cell_types_wReps]
sc_samples_wReps_cell_types <- do.call(rbind, strsplit(sc_samples_wReps, '__', 2))[,1]

bulk_data_use <- bulk_data[,samples_wReps]
CeNGEN_aggr_bio_reps_use <- CeNGEN_aggr_bio_reps[,sc_samples_wReps]


## subset down to just the samples with > 1 replicate in bulk and single cell

combined_cell_types <- intersect(unique(samples_wReps_cell_types), unique(sc_samples_wReps_cell_types))

bulk_data_use <- bulk_data_use[, str_split_fixed(colnames(bulk_data_use), 'r', 2)[,1] %in% combined_cell_types]
CeNGEN_aggr_bio_reps_use <- CeNGEN_aggr_bio_reps_use[, str_split_fixed(colnames(CeNGEN_aggr_bio_reps_use), '__', 2)[,1] %in% combined_cell_types]



## normalize bulk data to gene length

bulk_data_use <- (bulk_data_use/bulk_meta[rownames(bulk_data_use), 'Length']) * 1000



### integrate them, and if desired, save them

many_integrations <- pblapply(seq(101,150,1), function(v){
  set.seed(v)
  test_int <- integrate_geometricMean_biorep(bulk_data_use, CeNGEN_aggr_bio_reps_use,
                          pseudocount = 0.1, bulk_sep = 'r', single_cell_sep = '__')
  return(test_int) })
print('integrated samples')

many_integrations[[1]][1:10,1:10]

for(g in seq(1,3,1)){
  write.csv(many_integrations1[[g]], paste0('integrated_matrix_', g, '.csv'))
}


