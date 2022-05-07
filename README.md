# charcoal

Implementation of the changepoint localization method via the complementary sketching algorithm for high-dimensional regression coefficients. 

# Description of files

## R package

In `./R/` and `./man/` folders. Can be installed via `devtools::install_github('gaofengnan/charcoal')` in `R`.


<!-----------------------
## Python code

In the `./python/` folder
* `compsket.py`: file for the main algorithms
* `realdata.py`: file for implementing the real data example in Gao and Wang (2020)
* `example.ipynb`: IPython Notebook for the real data example

## MATLAB code (with a possible parallel computing implementation)

In the `./matlab/` folder
* `complementarySketching.m`: function for the main testing algorithm
* `differentialNetworkAnalysis.m`: specialized function for the nodewise regression testing on the gene interaction network example
* `main.m`: the script file processing the attached dataset
* `CD4_goodTREG_in_thymus.mat`: the preprocessed data for Matlab as in `./data/` 
## Data

* `CD4_TREG_in_thymus.csv`: preprocessed data for the real data example in Section 5 of Gao and Wang (2020). 
-->

# Reference
Gao, F. and Wang, T. (2022) Sparse changepoint detection in high-dimensional linear regression. Work-in-progress.
