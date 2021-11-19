## ncRNA
library(wbData)
library(edgeR)
library(MatrixGenerics)
library(reshape)
library(ComplexHeatmap)
library(dendsort)
library(tidyverse)
library(ggbeeswarm)
### PEM function ----


#adapted from Kryuchkova-Mostacci, et al., 207

# input is a genes x cell types matrix. I used the arithmetic mean across samples within cell types
fPem <- function(x){
  if(!all(is.na(x)))
  {
    x <- as.matrix(x)
    x[x<0] <- NA
    x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
    x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
    x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score
    
    x[x<1] <- 1
    x <- log10(x)
    
    x<- abs(x)				
    res <- x[-nrow(x),-ncol(x)]
    res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
  } else {
    res <- NA
    print("No data avalable.")
  }
  return(res)
}


### load data, normalize ----
bulk_data <- read.table('~/Bioinformatics/bsn5/bsn5_counts_110121.tsv')
bulk_meta <- read.table('~/Bioinformatics/bsn5/bsn5_bulk_counts_metadata.tsv')

NNLS_reduce <- read.table('NNLS_average_across_100_bootstraps.30Cells.111521.tsv')


bulk_data$ADFr99 <- NULL
bulk_data$M4r117 <- NULL

bulk_meta <- read.table('~/Bioinformatics/bsn5/bsn5_bulk_counts_metadata.tsv')

bulk_meta <- bulk_meta[rownames(bulk_data),]

bulk_data_pk <- (bulk_data/bulk_meta$Length) * 1000

bulk_raw_GeTMM <- DGEList(bulk_data_pk, group = str_split_fixed(colnames(bulk_data_pk), 'r', 2)[,1])
bulk_raw_GeTMM <- calcNormFactors(bulk_raw_GeTMM)
bulk_raw_GeTMM <- cpm(bulk_raw_GeTMM, normalized.lib.sizes = T)


### get average profile per cell type ----

aggr_raw_GeTMM <- bulk_raw_GeTMM
colnames(aggr_raw_GeTMM) <-str_split_fixed(colnames(aggr_raw_GeTMM),"r",2)[,1]
aggr_raw_GeTMM <- data.frame(vapply(unique(colnames(aggr_raw_GeTMM)), function(x) 
  rowMeans(aggr_raw_GeTMM[,colnames(aggr_raw_GeTMM)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(aggr_raw_GeTMM)) ))


## get reference, subset to ncRNAs of interest ----


ws277 <- wb_load_gene_ids('277') %>% data.frame(row.names = .$gene_id, .)

unique(ws277$biotype)

ws277_allnc <- ws277[ws277$biotype != 'gene' & ws277$biotype != 'protein_coding_gene' & 
                       ws277$biotype != 'miRNA_gene' & ws277$biotype != 'piRNA_gene' &
                       ws277$biotype != 'antisense_lncRNA_gene' & ws277$biotype != 'transposable_element_gene' &
                       ws277$biotype != 'rRNA_gene' & ws277$biotype != 'scRNA_gene', ]

aggr_nc_raw_GeTMM <- aggr_raw_GeTMM[rownames(aggr_raw_GeTMM) %in% rownames(ws277_allnc) ,]


#### load in neuron annotation metadata ----

neuron_meta.df <- read.table('../Neuron_metadata_hammarlund.tsv', sep = '\t',
                             header = T, row.names = 1)
neuron_meta.df$Modality_collapsed <- ' '
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Interneuron==1] = 'Interneuron'
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Unknown==1] = 'Unknown'
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Sensory==1] = 'Sensory'
neuron_meta.df$Modality_collapsed[neuron_meta.df$Modality..Motor==1] = 'Motor'



neuron_meta.df$Polymodal <- 'non'
neuron_meta.df$Polymodal[
  rowSums(neuron_meta.df[,c("Modality..Sensory", "Modality..Motor", "Modality..Interneuron", "Modality..Unknown")]) > 1] <- 'Polymodal'


neuron_meta_cut.df <- neuron_meta.df[colnames(aggr_nc_raw_GeTMM),]

### calculate PEM scores, genes x cell types

GeTMM_PEM <- fPem(aggr_raw_GeTMM)
GeTMM_PEM <- data.frame(GeTMM_PEM)

GeTMM_PEM_ncRNA <- GeTMM_PEM[intersect(rownames(GeTMM_PEM), rownames(ws277_allnc)),]

### calculate correlation to contaminants

GeTMM_cor_s_contaminant <- t(pbsapply(rownames(bulk_raw_GeTMM), function(gene){
  apply(NNLS_reduce[2:nrow(NNLS_reduce),], 1, function(tissue){
    stats::cor.test(bulk_raw_GeTMM[gene,], as.numeric(tissue), method = 'spearman')$estimate
  })
}))


GeTMM_cor_s_contaminant <- na.omit(GeTMM_cor_s_contaminant)



max_GeTMM_contaminant_corr <- data.frame(row.names = rownames(GeTMM_cor_s_contaminant),
                                         MatrixGenerics::rowMaxs(as.matrix(GeTMM_cor_s_contaminant)))



max_GeTMM_contaminant_corr_noncoding_genes <- intersect(rownames(max_GeTMM_contaminant_corr),
                                                        rownames(ws277_allnc))
max_GeTMM_contaminant_corr_noncoding <- data.frame(row.names = max_GeTMM_contaminant_corr_noncoding_genes,
                                                   ncRNA_corr = max_GeTMM_contaminant_corr[max_GeTMM_contaminant_corr_noncoding_genes,1]) 
max_GeTMM_contaminant_corr_noncoding


#### calculate binary expression using static threshold ----

nc_GeTMM_bin <- sapply(unique(str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]), function(cell){
  bulk <- bulk_raw_GeTMM[,str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]==cell]
  bulk <- bulk[intersect(rownames(bulk),
                         rownames(ws277_allnc)),]
  bin <- bulk > 5
  bin <- Matrix::rowMeans(bin)
  return((bin > 0.65) * 1)
})

ncRNA_expressed_somewhere <- rownames(nc_GeTMM_bin[Matrix::rowSums(nc_GeTMM_bin) > 0,])

#### maximum correlation graph for ncRNAs that are called expressed somemwhere ----


data.frame(ncRNA_corr = max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,]) |>
  ggplot() + geom_density( aes(x = ncRNA_corr, y = ..density..), color = '#B399D4', fill = '#B399D4', alpha = 0.9) +
  xlab('ncRNA Max correlation to contaminant') + geom_vline(xintercept = 0.2, size = 2, color = 'darkred') + 
  theme_classic(base_size = 20, base_line_size = 1) ### this one

length(max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,])
sum(max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,] > 0.2)
sum(max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,] > 0.2)/
  length(max_GeTMM_contaminant_corr_noncoding[ncRNA_expressed_somewhere,])



genes_low_max_contaminant_corr <- rownames(max_GeTMM_contaminant_corr)[max_GeTMM_contaminant_corr[,1] < 0.2]



#### genes passing thresholds ----
ncRNA_corr_cut.df <- pbsapply(colnames(aggr_nc_raw_GeTMM), function(cell){
  
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]
  
  bulk <- bulk_raw_GeTMM[, samples_cell_type %in% c(cell)]
  
  
  sapply(unlist(unique(ws277_allnc$biotype)), function(ncRNA){
    #print(genes)
    genes <- rownames(ws277_allnc[ws277_allnc$biotype==ncRNA,])
    genes <- intersect(genes, rownames(aggr_nc_raw_GeTMM))
    genes <- intersect(genes, genes_low_max_contaminant_corr)
    #print('genes')
    bulk_bin <- bulk[genes,] > 5
    
    bulk_bin_ave <- rowMeans(bulk_bin)
    
    sum(bulk_bin_ave > 0.65)
    
    
  })
})
ncRNA_corr_cut.df
reshape::melt(ncRNA_corr_cut.df) %>% #### this one
  data.frame(ncRNA_type = .[,1],
             cell_type = .[,2],
             total = .[,3]) %>%
  ggplot() + geom_col(aes(x = cell_type, y = total, fill = ncRNA_type)) + 
  theme_classic(base_size = 20) + xlab('Cell Type') +
  theme(axis.text.x = element_text(face = 'bold', angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(face = 'bold', color = 'black'),
        axis.title = element_text(face = 'bold', color = 'black'))


ncRNA_corr_cut.list.df <- lapply(colnames(aggr_nc_raw_GeTMM), function(cell){
  samples_cell_type <- str_split_fixed(colnames(bulk_raw_GeTMM), 'r', 2)[,1]
  bulk <- bulk_raw_GeTMM[, samples_cell_type %in% c(cell)]
  
  genes <- intersect(rownames(aggr_nc_raw_GeTMM), genes_low_max_contaminant_corr)
  bulk <- bulk[genes,]
  
  bulk_bin <- bulk[genes,] > 5
  
  bulk_bin_ave <- rowMeans(bulk_bin)
  
  
  return(names(bulk_bin_ave[bulk_bin_ave > 0.65]))
})
names(ncRNA_corr_cut.list.df) <- colnames(aggr_nc_raw_GeTMM)
ncRNA_corr_cut.list.df

ncRNA_corr_cut.list.df_table <- table(unlist(ncRNA_corr_cut.list.df))
ncRNA_corr_cut.list.df_table[ncRNA_corr_cut.list.df_table > 40]

sum(ncRNA_corr_cut.list.df_table >= 37)


data.frame(hist = seq(1,41,1),
           totals = sapply(seq(1,41,1), function(threshold){
  sum(ncRNA_corr_cut.list.df_table == threshold)
})) %>% ggplot() + geom_col(aes(x = hist, y = log10(totals)), color = 'black', fill = 'grey') + 
  scale_y_continuous(breaks = seq(0,3,1), labels = c('0', '10', '100', '1000')) +
  geom_vline(xintercept = 36.5, col = 'darkred', size = 1) +
  xlab('Total Cell types expressing ncRNA') + ylab('total ncRNAs') +
  theme(panel.background = element_blank(), axis.line = element_line(),
        axis.text = element_text(face = 'bold', size = 20), axis.title = element_text(face = 'bold', size = 30)
        )

ncRNA_corr_cut.list.df_table
highly_specific_ncRNAS <- aggr_nc_raw_GeTMM[names(ncRNA_corr_cut.list.df_table),][apply(GeTMM_PEM_ncRNA[names(ncRNA_corr_cut.list.df_table),], 1, max) > 0.65,]
highly_specific_ncRNAS <- highly_specific_ncRNAS[rowSums(highly_specific_ncRNAS > 3) <= 10,]
pan_ncRNAs <- aggr_nc_raw_GeTMM[names(ncRNA_corr_cut.list.df_table[ncRNA_corr_cut.list.df_table >=37]),]

dim(highly_specific_ncRNAS)

#pan_ncRNAs

library(circlize)


## set up colors ----
col_fun = colorRamp2(c(0, 3), c("#DEDCDC", "#Bb361d"))

col_fun2 = colorRamp2(c(0, 5), c("#DEDCDC", "#Bb361d"))

Heatmap(log10(1+pan_ncRNAs), show_row_names = F, name = ' ',
        col = col_fun2, show_row_dend = F, show_column_dend = F, column_order = column_order(hs_ncRNA_hmap))




  
ann <- data.frame(neuron_meta_cut.df$Modality_collapsed, neuron_meta_cut.df$Polymodal)
colnames(ann) <- c('Modality', 'Polymodal')
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = list(Modality = c(Interneuron = '#237400', 
                                                    Motor = '#D8AE00', 
                                                    Sensory = '#EE8366', 
                                                    Unknown = '#FCCCD8'),
                                       Polymodal = c(non = 'grey',
                                                     Polymodal = 'blue')),
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

dim(highly_specific_ncRNAS)

tmp.df <- t(scale(t(log10(1+highly_specific_ncRNAS))))
row_dend = dendsort(hclust(dist(tmp.df)))
col_fun = colorRamp2(c(0, 3), c("#DEDCDC", "#Bb361d"))
#col_fun = colorRamp2(c(min(unlist(tmp.df)),0, max(unlist(tmp.df))), c('blue', "#DEDCDC", "#Bb361d"))



### plot the highly specific genes ----
hmap <- Heatmap(
  log10(1+highly_specific_ncRNAS),
  #tmp.df,
  name = " ",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = T,
  cluster_rows = T,
  cluster_columns = F,
  column_order = order(neuron_meta_cut.df$Modality_collapsed),
  show_column_dend = F,
  show_row_dend = F,
  column_names_gp = gpar(fontsize = 15, fontface = 'bold'),
  top_annotation=colAnn)

draw(hmap, heatmap_legend_side="right", annotation_legend_side="right")


### plot the "pan-neuronal" genes ----
hmap2 <- Heatmap(
  log10(1+pan_ncRNAs),
  col = col_fun2,
  name = " ",
  show_row_names = FALSE,
  show_column_names = T,
  cluster_rows = T,
  cluster_columns = F,
  column_order = order(neuron_meta_cut.df$Modality_collapsed),
  show_column_dend = F,
  show_row_dend = F,
  column_names_gp = gpar(fontsize = 15, fontface = 'bold'),
  top_annotation=colAnn)

draw(hmap2, heatmap_legend_side="right", annotation_legend_side="right")




nrow(highly_specific_ncRNAS)


### plot max PEM specificity scores ----
apply(GeTMM_PEM_ncRNA[rowSums(GeTMM_PEM_ncRNA > 3) < 10 & rownames(GeTMM_PEM_ncRNA) %in% names(ncRNA_corr_cut.list.df_table),],
      1, max) %>% data.frame(max_PEM = .) %>% 
  ggplot() + geom_histogram(aes(x = (max_PEM)), color = 'black', fill = 'grey', binwidth = 0.005) +
  geom_vline(xintercept = (0.65), size = 2, color = 'darkred') +
  theme_classic(base_size = 20) +
  xlab('Max PEM specificity score')



#### how many highly specific genes per cell type


Specific_ncRNA_per_celltype <- data.frame(ncRNA = sapply(colnames(highly_specific_ncRNAS), function(cell){
  specific_genes <- rownames(highly_specific_ncRNAS)
  sum(GeTMM_PEM_ncRNA[specific_genes,cell] > 0.6)
  }),
           cell = names(highly_specific_ncRNAS_expressed),
           modality = sapply(names(highly_specific_ncRNAS_expressed), function(cell){
             neuron_meta_cut.df[cell, 'Modality_collapsed']
           })) 


#### plot the specific cells per cell type, grouped by modality ----
ggplot(Specific_ncRNA_per_celltype) + 
  geom_boxplot(aes(x = modality, y = (ncRNA))) + 
  geom_point(aes(x = modality, y = (ncRNA))) +
  theme_classic(base_size = 20) +
  xlab('') + ylab('Specific ncRNAs per cell type')


Sensory <- Specific_ncRNA_per_celltype |> filter(modality =='Sensory')
Sensory <- Sensory[,'ncRNA']
Interneuron <- Specific_ncRNA_per_celltype |> filter(modality =='Interneuron')
Interneuron <- Interneuron[,'ncRNA']
Motor <- Specific_ncRNA_per_celltype |> filter(modality =='Motor')
Motor <- Motor[,'ncRNA']


kruskal.test(x = c(Sensory, Interneuron), g = c(rep('Sensory', length(Sensory)), rep('Interneuron', length(Interneuron))))
kruskal.test(x = c(Motor, Interneuron), g = c(rep('Motor', length(Motor)), rep('Interneuron', length(Interneuron))))
kruskal.test(x = c(Sensory, Motor), g = c(rep('Sensory', length(Sensory)), rep('Motor', length(Motor))))





#### plot the number of pan neuronal and highly specific genes by RNA classs, pie charts ----

RNA_Classes <- unique(ws277_allnc$biotype)
specific_genes <- rownames(highly_specific_ncRNAS)


specific_genes_class_proportions <- sapply(RNA_Classes, function(class){
  sum(ws277_allnc[specific_genes, 'biotype'] == class)/length(specific_genes)
})

specific_genes_class_proportions <- data.frame(proportions = specific_genes_class_proportions,
           group = names(specific_genes_class_proportions))

pct <- round(specific_genes_class_proportions$proportions*100, digits = 1)
lbls <- specific_genes_class_proportions$group
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(specific_genes_class_proportions$proportions, labels = lbls, col=rainbow(length(lbls)),
    main="RNA classes in Highly Specific noncoding RNAs")


pan_genes <- rownames(pan_ncRNAs)
pan_genes_class_proportions <- sapply(RNA_Classes, function(class){
  sum(ws277_allnc[pan_genes, 'biotype'] == class)/length(pan_genes)
})

pan_genes_class_proportions <- data.frame(proportions = pan_genes_class_proportions,
                                               group = names(pan_genes_class_proportions))


pct <- round(pan_genes_class_proportions$proportions*100, digits = 1)
lbls <- pan_genes_class_proportions$group
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(pan_genes_class_proportions$proportions, labels = lbls, col=rainbow(length(lbls)),
    main="RNA classes in Pan-Neuronal noncoding RNAs")






