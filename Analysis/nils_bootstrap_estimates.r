#### Make 50 subsampled nnls estimates for the bulk samples ----

## note: use setseed for consistency


library(Seurat)
library(pbapply)
library(tidyverse)
library(dplyr)
library(nnls)
library(stringr)
library(ComplexHeatmap)


### seurat object downloaded from CeNGEN website downloads page, ~0.5 gb disk space, ~ 2gb RAM to open
sc_object <- readRDS('~/Bioinformatics/single_cell_data/100720_L4_all_cells_Seurat.rds')


sc_object <- sc_object[,sc_object$Tissue != 'Unknown' & sc_object$Tissue != 'Unannotated']


sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DB01')] <- 'DB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('DA9')] <- 'DA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VC_4_5')] <- 'VC'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VB01', 'VB02')] <- 'VB'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RMD_DV', 'RMD_LR')] <- 'RMD'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('RME_DV', 'RME_LR')] <- 'RME'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('VA12')] <- 'VA'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('IL2_DV', 'IL2_LR')] <- 'IL2'
sc_object@meta.data$Cell.type[sc_object@meta.data$Cell.type %in% c('AWC_ON', 'AWC_OFF')] <- 'AWC'
sc_object <- sc_object[,!(sc_object$Cell.type %in% c('RIV_stressed', 'SMD_stressed'))]

non_neuronal_list <- c('Excretory', 'Glia', 'Hypodermis', 'Intestine', 'Muscle_mesoderm', 'Pharynx', 'Reproductive')

sc_object$neuron_level <- sc_object$Tissue

sc_object$neuron_level[sc_object$neuron_level=='Neuron'] <- sc_object$Cell.type[sc_object$neuron_level=='Neuron']


sc_size <- sapply(unique(sc_object$neuron_level), function(cell){
  sum(sc_object$neuron_level == cell)
})

sc_size <- sc_size[order(sc_size)]
sc_size
## subset to 12 cells maximum


### cut down the size of the Seurat object to dramatically speed up the for loop
sc_object_cut <- CreateSeuratObject(counts = sc_object@assays$RNA@counts, 
                                    meta.data = sc_object@meta.data[,c('neuron_level', 'Cell.type', 'orig.ident', 'Barcode')])

NNLS_30_list <- pblapply(seq(101,200,1), function(seed){
  set.seed(seed)
  to_keep <- lapply(unique(sc_object_cut$neuron_level), function(cell){ ### list of cell barcodes to keep
    counter <- sum(sc_object_cut$neuron_level==cell)
    #print(cell)
    if(counter >= 30){
      keep_list <- sample(x = sc_object_cut$Barcode[sc_object_cut$neuron_level==cell], size = 20, replace = F)
      return(keep_list)
    }
    else { keep_list <- sc_object_cut$Barcode[sc_object_cut$neuron_level==cell]
    return(keep_list)} }) %>% unlist(.)
  
  res <- sc_object_cut[,to_keep]
  
  Idents(object = res) <- 'neuron_level'
  CeNGEN_max_arithMean <- sapply(unique(res$neuron_level), function(cell){
    #print(cell)
    #print(dim(res@assays$RNA@counts[,res$neuron_level==cell]))
    temp_mat <- res@assays$RNA@counts[,res$neuron_level==cell]
    temp_mat <- sweep(temp_mat, 2, colSums(temp_mat), '/')
    return(Matrix::rowMeans(temp_mat))
  })
  CeNGEN_max_arithMean <- CeNGEN_max_arithMean * 1000000
  
  
  CeNGEN_max_arithMean <- CeNGEN_max_arithMean[Matrix::rowSums(CeNGEN_max_arithMean > 10) > 0,]
  
  CeNGEN_max_arithMean_cv <- apply(CeNGEN_max_arithMean, 1, sd)
  CeNGEN_max_arithMean_cv <- CeNGEN_max_arithMean_cv/Matrix::rowMeans(CeNGEN_max_arithMean)
  
  
  CeNGEN_max_arithMean_cut <- CeNGEN_max_arithMean[CeNGEN_max_arithMean_cv > 3,]
  
  #sum(CeNGEN_12max_arithMean_cv > 10, na.rm = T)
  CeNGEN_max_arithMean_CPM.df <- data.frame(CeNGEN_max_arithMean_cut)
  CeNGEN_max_arithMean_CPM.df <- CeNGEN_max_arithMean_CPM.df[order(rownames(CeNGEN_max_arithMean_CPM.df)), 
                                                             order(colnames(CeNGEN_max_arithMean_CPM.df))]
  
  CeNGEN_max_arithMean_CPM.df$DD <- CeNGEN_max_arithMean_CPM.df$VD_DD
  CeNGEN_max_arithMean_CPM.df$VD <- CeNGEN_max_arithMean_CPM.df$VD_DD
  
  common.genes <- intersect(rownames(CeNGEN_max_arithMean_CPM.df), rownames(bulk_data))
  
  est <- sapply(colnames(bulk_data), function(sample1){
    cell <- str_split_fixed(sample1, 'r', 2)[,1]
    cell_types <- c(cell, non_neuronal_list)
    
    CeNGEN_CPMcut <- CeNGEN_max_arithMean_CPM.df[common.genes, cell_types]
    
    
    
    
    bulk <- as.numeric(bulk_data[common.genes, sample1])
    names(bulk) <- common.genes
    
    cell_types
    
    nnls(b=sqrt(as.matrix(bulk)), A=sqrt(as.matrix(CeNGEN_CPMcut)))$x
    
  })
  
  rownames(est) <- c('Neuron', non_neuronal_list)
  est <- sweep(est, 2, colSums(est), '/')
  
  return(est)
  
})


NNLS_reduce <- Reduce('+', NNLS_30_list)/length(NNLS_30_list)

write.table(NNLS_reduce, 'NNLS_average_across_100_bootstraps.30Cells.111521.tsv')


col_fun <- circlize::colorRamp2(c(0,1), c('white', 'black'))
Heatmap(NNLS_reduce, cluster_rows = F, cluster_columns = F, col = col_fun)
Heatmap(full_NNLS, cluster_rows = F, cluster_columns = F, col = col_fun)


bulk_data <- read.table('bsn5_counts_110121.tsv')

sc_size <- sapply(unique(sc_object$Cell.type[sc_object$Tissue=='Neuron']), function(cell){
  return(sum(sc_object$Cell.type==cell))
})
sc_size['DD'] <- sc_size['VD_DD']
sc_size['VD'] <- sc_size['VD_DD']

sc_size <- sc_size[order(names(sc_size))]
sc_size <- sc_size[-length(sc_size)]


tester <- sapply(colnames(bulk_data), function(sample1){
  cell <- str_split_fixed(sample1, 'r', 2)[,1]
  cell_types <- c(cell, non_neuronal_list)
  genes <- intersect(common.genes, names(CeNGEN_12max_arithMean_cv[CeNGEN_12max_arithMean_cv > 1]))
  
  CeNGEN_CPMcut <- CeNGEN_12max_arithMean_CPM.df[common.genes, cell_types]
  
  
  
  
  bulk <- as.numeric(bulk_data[common.genes, sample1])
  names(bulk) <- common.genes
  
  cell_types
  
  nnls(b=sqrt(as.matrix(bulk)), A=sqrt(as.matrix(CeNGEN_CPMcut)))$x
  
})







data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]]),
           Neuron = unlist(NNLS_reduce['Neuron',])) %>%
  ggplot() + 
  theme_classic(base_size = 20) + 
  geom_point(aes(x = sc_size, y = Neuron), color = 'black', size = 4) +
  geom_smooth(aes(x = sc_size, y = Neuron), method = 'lm', se = F) +
  xlab('Single Cell Cluster Size, Log10') +  ylim(0,1)


sc_size[str_split_fixed(colnames(NNLS_reduce1), 'r', 2)[,1]]
plot(log10(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]]), NNLS_reduce1['Neuron',], ylim = c(0,1))
plot(log10(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]]), full_NNLS['Neuron',], ylim = c(0,1))

plot(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]], NNLS_reduce['Neuron',], ylim = c(0,1))

plot(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]], full_NNLS['Neuron',])

full_NNLS

non_neuronal_list <- c("Excretory", "Glia", "Hypodermis", "Intestine", "Muscle_mesoderm", "Pharynx", "Rectal_cells", "Reproductive")




tester_log <- sapply(colnames(bulk_data), function(sample1){
  cell <- str_split_fixed(sample1, 'r', 2)[,1]
  cell_types <- c(cell, non_neuronal_list)
  
  CeNGEN_histNorm_cut <- CeNGEN_500max_arithMean_CPM.df[common.genes, cell_types]
  
  bulk <- as.numeric(bulk_data[common.genes, sample1])
  names(bulk) <- common.genes
  
  cell_types
  
  nnls(b=log2(1+as.matrix(bulk)), A=log2(1+as.matrix(CeNGEN_histNorm_cut)))$x
  
})
rownames(tester_log) <- c('Neuron', non_neuronal_list)
tester_log <- sweep(tester_log, 2, colSums(tester_log), '/')


tester[,1:10]
tester_log[,1:10]






CeNGEN_TPM <- read.table('~/Bioinformatics/single_cell_data/CeNGEN_TPM_080421.tsv')

CeNGEN_TPM_cv <- apply(CeNGEN_TPM, 1, function(gene){
  return(sd(gene)/mean(gene))
})
sum(CeNGEN_TPM_cv > 1, na.rm = T)

common.genes <- intersect(rownames(CeNGEN_TPM), rownames(bulk_data))

full_NNLS <- pbsapply(colnames(bulk_data), function(sample1){
    cell <- str_split_fixed(sample1, 'r', 2)[,1]
    cell_types <- c(cell, non_neuronal_list)
    genes <- intersect(common.genes, names(CeNGEN_TPM_cv[CeNGEN_TPM_cv > 5]))
    
    CeNGEN_CPMcut <- CeNGEN_TPM[common.genes, cell_types]
    
    
    
    
    bulk <- as.numeric(bulk_data[common.genes, sample1])
    names(bulk) <- common.genes
    
    cell_types
    
    nnls(b=sqrt(as.matrix(bulk)), A=sqrt(as.matrix(CeNGEN_CPMcut)))$x
    })
  
full_NNLS

rownames(full_NNLS) <- c('Neuron', non_neuronal_list)
full_NNLS <- sweep(full_NNLS, 2, colSums(full_NNLS), '/')
  


data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]]),
           Neuron = unlist(full_NNLS['Neuron',])) %>%
  ggplot() + 
  theme_classic(base_size = 20) + 
  geom_point(aes(x = sc_size, y = Neuron), color = 'black', size = 4) +
  geom_smooth(aes(x = sc_size, y = Neuron), method = 'lm', se = F) + 
  xlab('Single Cell Cluster Size, Log10') +  ylim(0,1)


data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]]),
           Neuron = unlist(NNLS_reduce1['Neuron',])) %>%
  ggplot() + 
  theme_classic(base_size = 20) + 
  geom_point(aes(x = sc_size, y = Neuron), color = 'black', size = 4) +
  geom_smooth(aes(x = sc_size, y = Neuron), method = 'lm', se = F) + 
  xlab('Single Cell Cluster Size, Log10') +  ylim(0,1)
ylab('Neuronal Proportion estimate') + ggtitle('12 cell bootstrap subsampled NNLS estimates')

full_lm <- lm(sc_size ~ Neuron, data = data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]]),
                                       Neuron = unlist(full_NNLS['Neuron',])))

downsampled_lm <- lm(sc_size ~ Neuron, data = data.frame(sc_size = log10(sc_size[str_split_fixed(colnames(NNLS_reduce), 'r', 2)[,1]]),
                                                         Neuron = unlist(NNLS_reduce1['Neuron',])))


summary(full_lm)
summary(downsampled_lm)











