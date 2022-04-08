###~~ protein_coding_genes detected somewhere, but not in certain neurons, raw bulk counts ----


# find a list of genes that are expressed in single cell somewhere -> main_list

## per cell type, list genes expressed cell_expressed, then use that to get the setdiff(main_list, cell_expressed) -> non_expressed
### subset it to exclude the non-neuronal ONLY genes
## in the bulk samples, find how many genes are expressed above the overall 5% FPR rate (non_neuronal) and ~14% FDR rate for neuronal ground truth
## 


## packages ----
library(edgeR)
library(pbapply)
library(tidyverse)
library(wbData)
library(cengenDataSC)

ws281 <- data.frame(wbData::wb_load_gene_ids('281'))
rownames(ws281) <- ws281$gene_id

protein_coding_genes <- ws281[ws281$biotype=='protein_coding_gene','gene_id']

protein_coding_genes <- intersect(protein_coding_genes, rownames(bulk_raw_GeTMM))

#cengen_sc_2_bulk

sc_gene_thresholds <- readRDS('~/Bioinformatics/single_cell_data/211028_genes_categorized_by_pattern.rds')

names(sc_gene_thresholds)
length(sc_gene_thresholds[[1]])

nondetected_threshold_genes <- sc_gene_thresholds[['nondetected']]
nondetected_threshold_genes <- intersect(nondetected_threshold_genes, polyA_genes)
non_neuronal_genes <- sc_gene_thresholds[['nonneuronal']]





### For each cell type that we have bulk data for: ----
###     make a list of genes that are: 1) Not detected in that cell in single cell at the liberal threshold
###                                    2) Are not labeled "undetected" in all single cell tissues
###                                    3) Are not detected exclusively in non-neuronal single cell clusters

sc_unexpr_2_list <- lapply(colnames(aggr_raw_GeTMM), function(cell){
  expr <- rownames(cengen_sc_1_bulk[cengen_sc_1_bulk[,cell] > 0,])
  unexpr <- setdiff(rownames(aggr_raw_GeTMM), expr)
  unexpr_protein <- intersect(unexpr, protein_coding_genes)
  unexpr_but_detected <- setdiff(unexpr_protein, nondetected_threshold_genes)
  unexpr_but_not_non_neuronal <- setdiff(unexpr_but_detected, non_neuronal_genes)
  
  return(unexpr_but_not_non_neuronal)
})
names(sc_unexpr_2_list) <- colnames(aggr_raw_GeTMM)
sapply(sc_unexpr_2_list, length)


### get a list of the potential genes, no duplicates
all_potential_unexpr <- unique(unlist(sc_unexpr_2_list))
length(all_potential_unexpr)




## non_neuronal GT ----
likely_non_neuronal <- read.table('~/Bioinformatics/Ground_truth_genesets/simpleMine_non_neuronal_genes_101521.tsv', sep = '\t',header = T)
likely_non_neuronal <- likely_non_neuronal[likely_non_neuronal$Putative.match.=='Yes',]

likely_non_neuronal_gt <- data.frame(row.names = likely_non_neuronal$WormBase.Gene.ID, 
                                     matrix(0, ncol=ncol(aggr_raw_GeTMM), nrow=nrow(likely_non_neuronal)))
colnames(likely_non_neuronal_gt) <- colnames(aggr_raw_GeTMM)



### exclude genes by correlation with contaminants
GeTMM_cor_s_contaminant_prc <- GeTMM_cor_s_contaminant[intersect(rownames(GeTMM_cor_s_contaminant), protein_coding_genes),]
GeTMM_cor_s_contaminant_prc_low <- GeTMM_cor_s_contaminant_prc[apply(GeTMM_cor_s_contaminant_prc, 1, max) < 0.3,]
dim(GeTMM_cor_s_contaminant_prc_low)
data.frame(prc_contam_cor_max = apply(GeTMM_cor_s_contaminant_prc[rownames(GeTMM_cor_s_contaminant_prc) %in% 
                                                                    all_potential_unexpr,], 1, max)) %>%
  ggplot() + geom_density(aes(x = prc_contam_cor_max))




#### plot all protein coding genes' correlations, and plot the correlations of the likely non-neuronal genes Figure 4A----
rbind(data.frame(. = apply(GeTMM_cor_s_contaminant_prc, 1, max), group = 'Non-Neuronal Genes'), 
      apply(GeTMM_cor_s_contaminant_prc[intersect(nn_genes, rownames(GeTMM_cor_s_contaminant_prc)),], 1, max) %>% 
        data.frame(., group = 'all protein coding genes') ) %>%
  ggplot() + geom_density(aes(x = ., color = group)) +
  xlim(-0.5,1) +
  theme_classic(base_size = 20) 
  


##### see how many non-neuronal genes are excluded by filtering by correlation with contaminants
nn_genes <- intersect(rownames(likely_non_neuronal_gt), common.genes)
nn_genes_low <- intersect(nn_genes, all_potential_unexpr) 
### subset to just the genes we're considering (many of these genes are filtered out by excluding non-neuronal genes & correlated genes)
nn_genes_low2 <- intersect(nn_genes, rownames(GeTMM_cor_s_contaminant_prc_low))

length(nn_genes_low)
length(nn_genes_low2)







##### set the threshold for calling genes expressed using an FDR to match the liberal threshold, see ground truth analysis script ----

Liberal_FDR_level <- 36 ### mean liberal fdr = 0.2, this is reached at a GeTMM threshold of 36 for the raw dataset



## calling the genes ----




GT_thresholded_FDR_L <- sapply(colnames(aggr_raw_GeTMM), function(cell){
  x1 <- rownames(aggr_raw_GeTMM[sc_unexpr_2_list[[cell]],][aggr_raw_GeTMM[sc_unexpr_2_list[[cell]], cell] > 36,])
  x2 <- rownames(GeTMM_cor_s_contaminant_prc_low)
  x3 <- intersect(x1, x2)
  length(x3)
})

GT_thresholded_FDR_L.gene.list <- sapply(colnames(aggr_raw_GeTMM), function(cell){
  x1 <- rownames(aggr_raw_GeTMM[sc_unexpr_2_list[[cell]],][aggr_raw_GeTMM[sc_unexpr_2_list[[cell]], cell] > 36,])
  x2 <- rownames(GeTMM_cor_s_contaminant_prc_low)
  x3 <- intersect(x1, x2)
  return(x3)
})

names(GT_thresholded_FDR_L.gene.list) <- colnames(aggr_raw_GeTMM)


sapply(GT_thresholded_FDR_L.gene.list, length)


GT_thresholded_FDR_L.df <- data.frame(protein_genes = as.numeric(GT_thresholded_FDR_L),
                                         cells = names(GT_thresholded_FDR_L))



##### bar plot of "new" genes Figure 4B----
GT_thresholded_FDR_L.df %>% ggplot(data = ., aes(x = cells, y = protein_genes)) + geom_col() +
  theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold')) +
  xlab('Cell Type') + ylab('"New" Genes Detected in Bulk\n(Sometimes Detected in Single Cell)')





new_genes.vs.size.df <- data.frame(raw_protein_genes = GT_thresholded_FDR_L.df$protein_genes,
           sc_size = sc_size[GT_thresholded_FDR_L.df$cells])


new_genes_decay <- nls((raw_protein_genes) ~ m + (M-m)*exp(-sc_size/alpha), data = new_genes.vs.size.df,
                       list(M = 350, alpha = 100, m = 18))


new_genes.df <- cbind(new_genes.vs.size.df, fit=predict(new_genes_decay))


ggplot(new_genes.df) +
  theme_classic(base_size = 25) +
  geom_line(aes(x=(sc_size), y=fit), color = "red3", linetype = "dashed") +
  geom_point(aes(x=(sc_size), y=(raw_protein_genes)), size = 4) +
  geom_hline(aes(yintercept = coef(new_genes_decay)[["m"]]), color = "gray50", linetype = "dotted") +
  #geom_smooth(aes(x=sc_size, y=raw_protein_genes), method = 'loess', se = F) + 
  xlab("Single Cell Cluster size") +
  ylab("'new' protein coding genes found in bulk")

names(GT_thresholded_FDR_0.15.gene.list) <- colnames(aggr_raw_GeTMM)



