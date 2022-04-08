# CeNGEN_integrated_analysis_biorxiv_code
code used analyzing data for the April 2022 CeNGEN biorxiv paper

Counts files with associated metadata are available here and on the cengen website (www.cengen.org)

The CeNGEN single cell seurat file is required for NNLS boostrapping.

Differential expression results files were not included for size limits, but can be recreated by running the files in the "differential expression" subdirectory in the "analysis" folder. The integrated analysis takes > 30 hours if run in one job, as it requires running the edgeR analysis 50 separate times.

Figure 1B: PCA_plots.R

Supp Figure 1B: single cell bulk cell type correlation.R

Figure 2B: TPR_FPR_FDR_calculation.R

Supp Figure 2: Normalization effects on non-neuronal FPR.R

Figure 3A-B: Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R then Accuracy_graphs.R

Figure 3C-D: Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R then MCC_graphs.R

Supp Figure 3C: Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R then Recall_TPR_graphs.R

Supp Figure 3D: Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R then Neuron_specificity_graphs.R

Supp Figure 3E: Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R then Accuracy_graphs.R

Supp Figure 3F: Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R then MCC_graphs.R

Supp Figure 3G: Loading_edgeR_results_and_calculating_TPR_FPR_FDR.R then nonNeuronal_Specificity_graphs.R

Figure 4A-C: Protein_Coding_Genes_Sometimes_Detected_in_SingleCell.R

Figure 4E: Protein_Coding_Genes_Never_Detected_in_SingleCell.R

Supp Figure 4A-B: NNLS_bootstrap_estimates.r

Supp Figure 4C-D: Protein_Coding_Genes_Sometimes_Detected_in_SingleCell

Figure 5 & Supp Figure 5: ncRNA_calling.r
