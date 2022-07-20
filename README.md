# RefTM
R package for RefTM
# Install
```
install.packages("devtools")
devtools::install_github("cuhklinlab/RefTM")
```
# Functions
This package includes following main functions:
- `RefTM` runs RefTM for the analysis of single-cell chromatin accessibility sequencing data. 
- `RefTM.LDA` runs RefTM-LDA with one specific number of topics.
- `RefTM.STM` runs RefTM-STM with one specific number of topics.
- `RefTM_postprocess` postprocess the result to obtain a final cell-by-topic matrix as the output of RefTM.
- `RefTM_tsne` visualize the output of RefTM via a t-SNE plot.
- `RefTM_motif` runs motif enrichment based on the output of RefTM.

# Documentation
Please check the [vigenette](https://github.com/cuhklinlab/RefTM/wiki) for a tutorial.
