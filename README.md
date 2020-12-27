# EMeth
 This package provides implementations of the EMeth algorithm. It contains two functions.
- _emeth_, which provides implementation of two families of EMeth, normal and laplace, accordingto what likelihood is used. Generally we recommend the use of laplace family.
- _cv.emeth_, which helps tuning the ridge penalty by cross validation.

## Prerequisite
 EMeth package requires R version 3.6.0 or higher. It requires the package ```quadprog```

## Installation 
 To install this package in R, use 
 
 ```
    library("devtools");
    install_github("Hanyuz1996/EMeth")
 ```
 
## Example
  An example use:
  ```
     library("EMeth")
     Emeth_laplace = cv.emeth(Y, eta, mu, aber = TRUE, V='c', init = 'default',
                               family = 'laplace', nu = penalty, folds = 5, 
                               maxiter = 50, verbose = TRUE)
  ```
where ```Y``` is the DNA methylation data matrix with rows for CpG probes and columns for samples, ```eta``` is a vector of tumor purity, and ```mu``` is reference of cell type-specific DNA methylation, with each column for one cell type. This function ```cv.emeth``` automatically runs the cross validation procedure and EMeth fitting with laplace . Example data are provided in the example folder. 
  
## Support
  You can contact zhanghy1996@gmail.com
  
## Citation
  If you use the software, please cite our paper: Zhang et al. 2021. The pipelines for simulation studies and real data analysis in this paper are contained in [this repository](https://github.com/Sun-lab/dMeth).
  
## Authors
  Hanyu Zhang (University of Washington)
  
  Wei Sun (Fred Hutchinson Cancer Research Center)
                              
## Reference

Zhang et al. 2021, EMeth: An EM algorithm for cell type decomposition based on DNA methylation data.


 
 
