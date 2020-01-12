# EMeth
This is a package for conducting cell type deconvolution by EMeth algorithm.  This package include 4 functions:

- deconv_EM: Perform EM algorithm with normal model.
- deconv_CV:  Perform cross validation for penalized OriEM.
- deconv_laplace: Perform EM algorithm with Laplace distribution model.
- deconv_laplace_CV: Perform  cross validation for penalized LaplaceEM.

The folder in this repository includes:

- EMeth: The four functions in EMeth.
- pipelines: The code for the three studies in the paper.
- source: utility functions for the three studies.