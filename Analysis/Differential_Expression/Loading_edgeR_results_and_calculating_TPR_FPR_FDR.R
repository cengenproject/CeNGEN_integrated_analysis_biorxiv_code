#### setting up contrasts for ground truth


### load libraries ----
library(tidyverse)
library(dplyr)
library(reshape)

library(wbData)
library(bayestestR)
library(edgeR)
library(pbapply)

library(ggbeeswarm)
library(ComplexHeatmap)
library(circlize)
library(dendsort)



### functions ----
get_tpr <- function(expression, truth, threshold, na.rm = TRUE){
  # True Positive Rate, aka sensitivity, aka recall
  # TPR = TP/(TP+FN) = TP/P
  bin <- expression >= threshold
  bin[is.na(bin)] <- 0
  return(sum(bin * truth)/sum(truth))
}
get_fpr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Positive Rate
  # FPR = FP/(FP+TN) = FP/N
  bin <- expression >= threshold
  bin[is.na(bin)] <- 0
  return(sum(bin * (!truth))/sum(!(truth)))
}
get_fdr <- function(expression, truth, threshold, na.rm = TRUE){
  # False Discovery Rate
  # FDR = FP/(FP+TP) = 1 - PPV
  bin <- expression >= threshold
  bin[is.na(bin)] <- 0
  fdr <- sum(bin * (!truth))/(sum(bin*(!truth)) + sum(bin*truth))
  if(is.nan(fdr))
    fdr <- 0
  return(fdr)
}


### load data

bulk_data <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.tsv')

CeNGEN_TPM <- read.table('CeNGEN_TPM_080421.tsv')


### load ground truth dataset ----

bulk_gt <- read.csv('bulk_neuronal_ground_truth_042021.csv', row.names = 1)


nonNeuronal_GT <- read.csv('bulk_nonNeuronal_ground_truth_040122.csv', header = 1)
rownames(nonNeuronal_GT)


ws281 <- data.frame(wb_load_gene_ids('281'))
rownames(ws281) <- ws281$gene_id

## subset to common genes

common.genes <- intersect(rownames(bulk_data), rownames(CeNGEN_TPM))
common.genes <- intersect(common.genes, ws281[ws281$biotype=='protein_coding_gene', 'gene_id'])

nn_genes <- intersect(nn_genes, common.genes)


### using the bulk dataset, generate a list of possible cell-cell pairs

cell_types <- unique(str_split_fixed(colnames(bulk_load), 'r', 2)[,1])

contrast_list <- list()
for(c in cell_types){
  for(g in cell_types[cell_types!=c]){
    new_contrast <- paste0(c,'-(',g,')')
    contrast_list <- c(contrast_list, new_contrast)
  }
}
names(contrast_list) <- contrast_list

cell_types_cut <- cell_types
short_contrast_list <- list()
for(c in cell_types_cut){
  cell_types_cut <- cell_types_cut[cell_types_cut!=c]
  for(g in cell_types_cut){
    new_contrast <- paste0(c,'-(',g,')')
    short_contrast_list <- c(short_contrast_list, new_contrast)
  }
}
names(short_contrast_list) <- short_contrast_list
length(short_contrast_list)
length(contrast_list)

contrast_list_match <- list()
for(c in cell_types){
  for(g in cell_types[cell_types!=c]){
    new_contrast <- paste0(c,'-(',g,')')
    
    if(new_contrast %in% names(short_contrast_list)){
      matched_contrast <- new_contrast
      contrast_list_match <- c(contrast_list_match, matched_contrast)
    }
    else{
      matched_contrast <- paste0(g,'-(',c,')')
      contrast_list_match <- c(contrast_list_match, matched_contrast)
    }
    #contrast_list <- c(contrast_list, new_contrast)
  }
}
names(contrast_list_match) <- contrast_list

matched_contrast <- data.frame( contrast = names(contrast_list_match),
                                match = unlist(contrast_list_match))



##  for each contrast, grab those two columns from the ground truth, in the order specified by the contrast ----


contrast_gt_list <- sapply(matched_contrast$contrast, function(contrast){
  cell1 <- str_split_fixed(contrast, '-', 2)[,1]
  cell2 <- str_split_fixed(contrast, '-', 2)[,2]
  cell2 <- gsub(x = cell2, pattern = '\\(|\\)', replacement = '')
  
  gt_sub <- bulk_gt[,c(cell1, cell2)]
  gt_sub <- gt_sub[rowSums(gt_sub) <= 1,]
  gt_return <- gt_sub[,1]
  names(gt_return) <- rownames(gt_sub)
  
  return(gt_return)
  
})

names(contrast_gt_list) <- matched_contrast$contrast




### load the edgeR tables, examples provided here, but the files are too large to share on github ----

raw_qlfs <- readRDS('~/Bioinformatics/bsn9/raw_qlf_tables_bsn9_112821.rds') ### uncorrected bulk differential expression tables

many_harmonized <- readRDS('~/Bioinformatics/bsn9/many_harmonized_edgeR_113021.rds') #### integrated differential expression tables

dim(many_harmonized[[1]])
dim(many_integrations[[1]])


## For each cell-cell pair, in each direction, calculate the True Positive Rate, False Positive Rate, and the False Discovery Rate

auc_TPR_consensus_p.hmp_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- many_harmonized[[contrast]]
    
    testing_gt1 <- contrast_gt_list[[contrast]]
    
    
    testing_harmony_true <- testing_harmony[names(testing_gt1),3]
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1), 2] < 2] <- NA
    
    TPR <- get_tpr(testing_harmony_true, testing_gt1, threshold = 40)
    
    return(TPR)
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- many_harmonized[[match]]
      
      testing_gt1 <- contrast_gt_list[[contrast]]
      
      
      testing_harmony_true <- testing_harmony[names(testing_gt1),3]
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1), 2] > (-2)] <- NA
      
      
      TPR <- get_tpr(testing_harmony_true, testing_gt1, threshold = 40)
      
      return(TPR)
      
    }
    else{NA}
  }
})

auc_TPR_raw_PValue_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- raw_qlfs[[contrast]]
    
    testing_gt1 <- contrast_gt_list[[contrast]]
    
    
    testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] < 2] <- NA
    
    TPR <- get_tpr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
    
    return(TPR)
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- raw_qlfs[[match]]
      
      testing_gt1 <- contrast_gt_list[[contrast]]
      
      
      testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] > (-2)] <- NA
      
      TPR <- get_tpr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
      
      return(TPR)
      
    }
    else{NA}
  }
})
auc_TPR_consensus_p.hmp_directional <- na.omit(auc_TPR_consensus_p.hmp_directional)
auc_TPR_raw_PValue_directional <- na.omit(auc_TPR_raw_PValue_directional)

auc_FPR_consensus_p.hmp_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- many_harmonized[[contrast]]
    
    testing_gt1 <- contrast_gt_list[[contrast]]
    
    
    testing_harmony_true <- testing_harmony[names(testing_gt1),3]
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1), 2] < 2] <- NA
    
    FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = 40)
    
    return(FPR)
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- many_harmonized[[match]]
      
      testing_gt1 <- contrast_gt_list[[contrast]]
      
      
      testing_harmony_true <- testing_harmony[names(testing_gt1),3]
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1), 2] > (-2)] <- NA
      
      
      FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = 40)
      
      return(FPR)
      
    }
    else{NA}
  }
})

auc_FPR_raw_PValue_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- raw_qlfs[[contrast]]
    
    testing_gt1 <- contrast_gt_list[[contrast]]
    
    
    testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] < 2] <- NA
    
    FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
    
    return(FPR)
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- raw_qlfs[[match]]
      
      testing_gt1 <- contrast_gt_list[[contrast]]
      
      
      testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] > (-2)] <- NA
      
      FPR <- get_fpr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
      
      return(FPR)
      
    }
    else{NA}
  }
})
auc_FPR_consensus_p.hmp_directional <- na.omit(auc_FPR_consensus_p.hmp_directional)
auc_FPR_raw_PValue_directional <- na.omit(auc_FPR_raw_PValue_directional)

auc_FDR_consensus_p.hmp_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- many_harmonized[[contrast]]
    
    testing_gt1 <- contrast_gt_list[[contrast]]
    
    
    testing_harmony_true <- testing_harmony[names(testing_gt1),3]
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1), 2] < 2] <- NA
    
    FDR <- get_fdr(testing_harmony_true, testing_gt1, threshold = 40)
    
    return(FDR)
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- many_harmonized[[match]]
      
      testing_gt1 <- contrast_gt_list[[contrast]]
      
      
      testing_harmony_true <- testing_harmony[names(testing_gt1),3]
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1), 2] > (-2)] <- NA
      
      
      FDR <- get_fdr(testing_harmony_true, testing_gt1, threshold = 40)
      
      return(FDR)
      
    }
    else{NA}
  }
})

auc_FDR_raw_PValue_directional <- pbsapply(matched_contrast$contrast, function(contrast){
  
  match <- matched_contrast[matched_contrast$contrast == contrast, 'match']
  
  if(contrast %in% names(many_harmonized)){
    testing_harmony <- raw_qlfs[[contrast]]
    
    testing_gt1 <- contrast_gt_list[[contrast]]
    
    testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
    names(testing_harmony_true) <- names(testing_gt1)
    
    testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] < 2] <- NA
    
    FDR <- get_fdr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
    
    return(FDR)
    
    
  }
  else{
    if(match %in%  names(many_harmonized)){
      
      testing_harmony <- raw_qlfs[[match]]
      
      testing_gt1 <- contrast_gt_list[[contrast]]
      
      testing_harmony_true <- -log10(testing_harmony[names(testing_gt1),'PValue'])
      names(testing_harmony_true) <- names(testing_gt1)
      
      testing_harmony_true[testing_harmony[names(testing_gt1),'logFC'] > -2] <- NA
      
      FDR <- get_fdr(testing_harmony_true, testing_gt1, threshold = -log10(0.05))
      
      
      return(FDR)
      
    }
    else{NA}
  }
})
auc_FDR_consensus_p.hmp_directional <- na.omit(auc_FDR_consensus_p.hmp_directional)
auc_FDR_raw_PValue_directional <- na.omit(auc_FDR_raw_PValue_directional)

auc_FDR_raw_PValue_directional

### calculate total true positives for accuracy score

consensus_p.hmp_total_TP <- sapply(names(auc_TPR_consensus_p.hmp_directional), function(contrast){
  gt <- contrast_gt_list[[contrast]]
  total_true <- sum(gt==1)
  return(total_true * auc_TPR_consensus_p.hmp_directional[[contrast]])
})
consensus_p.hmp_total_TN <- sapply(names(auc_FPR_consensus_p.hmp_directional), function(contrast){
  gt <- contrast_gt_list[[contrast]]
  total_true <- sum(gt==0)
  return(total_true * (1-auc_FPR_consensus_p.hmp_directional[[contrast]]))
})


raw_PValue.hmp_total_TP <- sapply(names(auc_TPR_raw_PValue_directional), function(contrast){
  gt <- contrast_gt_list[[contrast]]
  total_true <- sum(gt==1)
  return(total_true * auc_TPR_raw_PValue_directional[[contrast]])
})
raw_PValue_total_TN <- sapply(names(auc_FPR_raw_PValue_directional), function(contrast){
  gt <- contrast_gt_list[[contrast]]
  total_true <- sum(gt==0)
  return(total_true * (1-auc_FPR_raw_PValue_directional[[contrast]]))
})

totals_per_contrast <- sapply(contrast_gt_list, length)

consensus_p.hmp_accuracy <- (consensus_p.hmp_total_TP + consensus_p.hmp_total_TN)/totals_per_contrast[names(consensus_p.hmp_total_TP)]
consensus_p.hmp_accuracy[is.na(consensus_p.hmp_accuracy)] <- 0

raw_PValue_accuracy <- (raw_PValue.hmp_total_TP + raw_PValue_total_TN)/totals_per_contrast[names(raw_PValue.hmp_total_TP)]
raw_PValue_accuracy[is.na(raw_PValue_accuracy)] <- 0

summary(raw_PValue_accuracy)
summary(consensus_p.hmp_accuracy)


### Matthew's Correlation Coefficient calculation

mcc <- function(tp, tn, gt){
  total_true <- sum(gt==1)
  total_false <- sum(gt==0)
  
  tp <- tp
  tn <- tn
  fn <- total_true - tp
  fp <- total_false - tn
  
  numer <- (tp * tn) - (fp * fn)
  denom <- sqrt((tp + fp)* (tp + fn) * (tn + fp) * (tn + fn)) 
  
  return(numer/denom)
}

raw_pValue_MCC <- sapply(names(auc_TPR_raw_PValue_directional), function(contrast){
  tp <- raw_PValue.hmp_total_TP[[contrast]]
  tn <- raw_PValue_total_TN[[contrast]]
  gt <- contrast_gt_list[[contrast]]
  mcc(tp, tn, gt)
})

consensus_MCC <- sapply(names(auc_TPR_consensus_p.hmp_directional), function(contrast){
  tp <- consensus_p.hmp_total_TP[[contrast]]
  tn <- consensus_p.hmp_total_TN[[contrast]]
  gt <- contrast_gt_list[[contrast]]
  mcc(tp, tn, gt)
})



summary(raw_pValue_MCC)
summary(consensus_MCC)



