# Chen

This is the pipelines for results in Chen data. Details of experimental design and resutls are provided in the supplementary materials. This pipeline includes:

- step0-GetData.R: get and preprocess the methylation data matrix.
- step1-estMeanVar.R: estimate the mean and variance from methylation and remove probes with extremely large variance. Then separate the data and obtain dat.gen to generate mixture and dat.est to construct the reference data.
- step2-filterProbes.R:  select probes based on dat.est and save in 'probes.RData'
- step3-RunMethy.R: generate methylation mixture and run cell type deconvolution. 


- step4-runExpr.R: generate DNA expression mixture data and then perform CIBERSORT and compare the result with results from methylation data.