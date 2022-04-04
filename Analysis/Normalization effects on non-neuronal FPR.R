#### comparing TMM and GeTMM normalized counts for non-neuronal FPR



CeNGEN_aggr_bio_reps <- read.table('~/Bioinformatics/single_cell_data/CeNGEN_aggr_bio_reps_101821.tsv')

many_integrations <- pblapply(seq(101,150,1), function(v){
  set.seed(v)
  test_int <- integrate_geometricMean_biorep(bulk_data, CeNGEN_aggr_bio_reps, pseudocount = 0.1, bulk_sep = 'r', single_cell_sep = '__')
  return(test_int)
})

bulk_raw_GeTMM <- DGEList(counts = bulk_data_pk)
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM, method = 'TMM')
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)

many_integrations_TMM <- pblapply(many_integrations, function(v){
  g <- DGEList(counts = v)
  g <- calcNormFactors(g, method = 'TMM')
  g <- cpm(g, normalized.lib.sizes = T)
  return(g)
})


bulk_raw_GeTMM <- DGEList(counts = bulk_data_pk)
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM, method = 'TMM')
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)

bulk_raw_TMM <- DGEList(counts = bulk_data)
bulk_raw_TMM <- calcNormFactors(bulk_raw_TMM, method = 'TMM')
bulk_raw_TMM <- cpm(bulk_raw_TMM, normalized.lib.sizes = T)

bulk_raw_GeTMM[1:5,1:10]
bulk_raw_TMM[1:5,1:10]


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


common.genes <- intersect(rownames(aggr_raw_TMM), rownames(CeNGEN_TPM))


neurons <- intersect(colnames(aggr_raw_TMM), colnames(all_genes_bulk))
genes <- intersect(rownames(aggr_raw_TMM), rownames(all_genes_bulk))
#genes <- intersect(genes, rownames(prop_by_type_adjusted))



aggr_raw_GeTMM_plot <- aggr_raw_GeTMM[genes, neurons]
aggr_raw_TMM_plot <- aggr_raw_TMM[genes, neurons]

#props_adjust_plot <- prop_by_type_adjusted[genes, neurons]

all_gt <- all_genes_bulk[genes, neurons]

###~~ calc TPR FPR FDR ----


diags_aggr_raw_TMM_plot <- tibble(threshold = c(0,2**seq(-4,12,0.1)),
                                              TPR = map_dbl(threshold, ~get_tpr(aggr_raw_TMM_plot, all_gt, .x)),
                                              FPR = map_dbl(threshold, ~get_fpr(aggr_raw_TMM_plot, all_gt, .x)),
                                              FDR = map_dbl(threshold, ~get_fdr(aggr_raw_TMM_plot, all_gt, .x)),
                                              counts = "aggr_raw_TMM_plot")


diags_aggr_raw_GeTMM_plot <- tibble(threshold = c(0,2**seq(-4,12,0.1)),
                                    TPR = map_dbl(threshold, ~get_tpr(aggr_raw_GeTMM_plot, all_gt, .x)),
                                    FPR = map_dbl(threshold, ~get_fpr(aggr_raw_GeTMM_plot, all_gt, .x)),
                                    FDR = map_dbl(threshold, ~get_fdr(aggr_raw_GeTMM_plot, all_gt, .x)),
                                    
                                    counts = "aggr_raw_GeTMM_plot")



bind_rows(diags_aggr_raw_TMM_plot,
          diags_aggr_raw_GeTMM_plot) %>%
  ggplot(aes(x = FPR, y=TPR, color= counts)) +
  geom_point(size = 3) +
  theme_classic(base_size = 20) 
  theme(legend.position = '') +
  xlim(0,0.99)




likely_non_neuronal <- read.table('~/Bioinformatics/Ground_truth_genesets/simpleMine_non_neuronal_genes_101521.tsv', sep = '\t',header = T)
likely_non_neuronal <- likely_non_neuronal[likely_non_neuronal$Putative.match.=='Yes',]

likely_non_neuronal_gt <- data.frame(row.names = likely_non_neuronal$WormBase.Gene.ID, 
                                     matrix(0, ncol=ncol(all_genes_bulk), nrow=nrow(likely_non_neuronal)))
colnames(likely_non_neuronal_gt) <- colnames(all_genes_bulk)

nn_genes <- intersect(rownames(likely_non_neuronal_gt), common.genes)


likely_non_neuronal_gt_cut <- likely_non_neuronal_gt[nn_genes, neurons]


aggr_raw_GeTMM_nnplot <- aggr_raw_GeTMM[nn_genes, neurons]
aggr_raw_TMM_nnplot <- aggr_raw_TMM[nn_genes, neurons]

dim(likely_non_neuronal_gt_cut)
dim(CeNGEN_TPM_nnplot)

#diags_props_adjust_nnplot <- tibble(threshold = seq(0,1,0.001),
#                                    FPR = map_dbl(threshold, ~get_fpr(prop_by_type_adjusted_nnplot, likely_non_neuronal_gt_cut, .x)),
#                                    counts = "prop_by_type_adjusted_nnplot")

diags_aggr_raw_TMM_nnplot <- tibble(threshold = c(0,2**seq(-3,10,0.1)),
                                                FPR = map_dbl(threshold, ~get_fpr(aggr_raw_TMM_nnplot, likely_non_neuronal_gt_cut, .x)),
                                                counts = "aggr_raw_TMM_nnplot")

diags_aggr_raw_GeTMM_nnplot <- tibble(threshold = c(0,2**seq(-3,10,0.1)),
                                      FPR = map_dbl(threshold, ~get_fpr(aggr_raw_GeTMM_nnplot, likely_non_neuronal_gt_cut, .x)),
                                      counts = "aggr_raw_GeTMM_nnplot")


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

data.frame(GeTMM_mean = apply(aggr_raw_TMM, 1, max),
           TMM_mean = apply(bulk_raw_GeTMM, 1, max)) %>% reshape::melt(.) %>%
  ggplot() +
  geom_density(aes(x = log10(value+1), fill = variable), alpha = 0.6) +
  
  scale_x_continuous(breaks = log10(c(1,11,101,1001, 10001)),
                     labels = c('0', '10', '100', '1,000', '10k'))+
  theme_classic(base_size = 35) +
  theme(legend.position = '') +
  xlab('Normalized Counts log10(+1)')



auc(log10(diags_aggr_raw_TMM_nnplot$threshold+1)/log10(max(diags_aggr_raw_TMM_nnplot$threshold+1)), diags_aggr_raw_TMM_nnplot$FPR)
auc(log10(diags_aggr_raw_GeTMM_nnplot$threshold+1)/log10(max(diags_aggr_raw_GeTMM_nnplot$threshold+1)), diags_aggr_raw_GeTMM_nnplot$FPR)





bulk_load
