# TCGA Real Data Analysis

This pipelines are codes to run real data analysis from TCGA data base. The files includes:

- COAD.R, LUAD.R, LUSC.R, SKCM.R: preprocess methylation data for each cancer type: generate a DNA methylation data matrix. Probes used are pre-stored in 'ref_966probes.RData', which are selected in purified cell type data in Combined Data Study. 
- Gene expression data are pre-estimated by cibersortx and saved in a different file.
- runCOAD.R, runLUAD.R, runLUSC.R, runSKCM.R: run EMeth and benchmark algorithms and compare with previously estimated result from CIBERSORTX. 

When running analysis, one can directly run the files starting with 'run*' for certain cancer type. 