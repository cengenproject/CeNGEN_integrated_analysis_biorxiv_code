# CeNGEN_integrated_analysis_biorxiv_code
code used analyzing data for the April 2022 CeNGEN biorxiv paper

Counts files with associated metadata are available here and on the cengen website (www.cengen.org)

The CeNGEN single cell seurat file is required for NNLS boostrapping.

Differential expression results files were not included for size limits, but can be recreated by running the files in the "differential expression" subdirectory in the "analysis" folder. The integrated analysis takes > 30 hours if run in one job, as it requires running the edgeR analysis 50 separate times.

