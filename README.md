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
     em1 = cv.emeth(Y, eta, mu, aber = TRUE, V='c', init = 'default',
     	family = 'laplace', nu = penalty, folds = 5, 
	maxiter = 50, verbose = TRUE)
     cell_type_prop_est = em1$result$rho
  ```
where 
* ```Y``` is the DNA methylation data matrix with rows for CpG probes and columns for samples, 
* ```eta``` is a vector of tumor purity that can be set as 0 for non-tumor studies, and 
* ```mu``` is the reference data for cell type-specific DNA methylation, with each column for one cell type. 
* ```nu``` is the penalty parameters, for example. we use ```nrow(Y)*(10^seq(-2,1,1)) ``` in our TCGA analysis. 

This function ```cv.emeth``` automatically runs the cross validation procedure. Please see the help document for this function ```?cv.emeth``` for more details of these and other parameters. Example data are provided in the example folder. 
 
The output is a list with three elements

* ```result```: The result of EMeth algorithm using the penalty value selected by cross-validation. It is a list itself. 
* ```choosenu```: The value of the ```nu``` (the penalty) chosen by the cross-validation.
* ```losslist```: A matrix saving the loss for each fold and each choice of ```nu```.

Most often the output of interest is ```result$rho```, which is a matrix of the cell type proportion estimates (rows for samples and columns for cell type). Another output that could be useful is ```result$gamma```, which is a matrix. The (i,j)-th entry of  ```result$gamma``` is the probability that the i-th probe in the j-th sampl is aberrant (i.e., the cell type-specific methylation is bulk sample is not consistent with those in reference).


## Support
  You can contact zhanghy1996@gmail.com
  
## Citation
  If you use the software, please cite our paper: Zhang et al. 2021. The pipelines for simulation studies and real data analysis in this paper are contained in [this repository](https://github.com/Sun-lab/dMeth).
  
## Authors
  Hanyu Zhang (University of Washington)
  
  Wei Sun (Fred Hutchinson Cancer Research Center)
                              
## Reference

Zhang et al. 2021, EMeth: An EM algorithm for cell type decomposition based on DNA methylation data.


 
 
