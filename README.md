# EMeth
 This package provides implementations of the EMeth algorithm. It contains two functions.
- _emeth_, which provides implementation of two families of EMeth, normal and laplace, accordingto what likelihood is used. Generally we recommend the use of laplace family.
- _cv.emeth_, which helps tuning the ridge penalty by cross validation.

## Prerequisite
 EMeth package requires R version 3.6.0 or higher. 

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
     Emeth_laplace = cv.emeth(Y,eta,mu,aber = TRUE, V='c', init = 'default',
                               family = 'laplace', nu = penalty, folds = 5, maxiter = 50, verbose = TRUE)
  ```
  which automatically runs the cross validation procedure and EMeth fitting with laplace . Example data are provided in the example folder. 
  
## Support
  You can contact zhanghy1996@gmail.com
  
## Citation
  If you use the software, please cite the following paper. _To be added: paper link_
  In this paper, the pipelines for simulation studies and real data analysis are contained in [this repository](https://github.com/Sun-lab/dMeth).
  
## Authors
  Hanyu Zhang (University of Washington)
  
  Wei Sun (Fred Hutchinson Cancer Research Center)
                              

 
 
