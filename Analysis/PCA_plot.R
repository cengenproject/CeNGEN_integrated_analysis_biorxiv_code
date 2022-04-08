### script for PCA plot, Figure 1B

## libraries ----

library(ggplot2)
library(DESeq2)


## load data & metadata ----
bulk_data <- read.table('Barrett_et_al_2022_CeNGEN_bulk_RNAseq_data.tsv')

neuron_meta.df <- read.table('~/Bioinformatics/Neuron_metadata_hammarlund.tsv', sep = '\t',
                             header = T, row.names = 1)
neuron_meta.df$Modality_collapsed <- ' '
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Interneuron==1] = 'Interneuron'
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Unknown==1] = 'Unknown'
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Sensory==1] = 'Sensory'
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Motor==1] = 'Motor'

## keep the seed reproducible
set.seed(12345)

bulk_data_deseqdata <- DESeq2::DESeqDataSetFromMatrix(bulk_data,
                                                      colData = data.frame(condition=str_split_fixed(colnames(bulk_data), 'r', 2)[,1]),
                                                      design = ~ condition)


bulk_data_deseqdata


## use DESeq based variance stabilization transformation
bulk_data_VST <- DESeq2::vst(bulk_data_deseqdata)

## calculate first two Principal Components, and add modality info to the resultant dataframe
plt <- DESeq2::plotPCA(bulk_data_VST, returnData = T)
plt$modality <- sapply(plt$group, function(cell){
  modality = neuron_meta.df[cell, 'Modality_collapsed']
})

## plot the PCA components with color coding by modality and labeling by cell type
ggplot() + geom_label(data = plt, aes(x = PC1, y = PC2, label = name, fill = modality, color = name), label.size = NA) +
  scale_color_manual(values = rep('white', 160)) + 
  theme_classic(base_size = 25) + theme(legend.position = '')
ggsave('PCA_plot.pdf')
