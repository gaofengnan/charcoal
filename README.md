# charcoal - A changepoint detection algorithm in high-dimensional linear regression

(By [Fengnan Gao](https://gaofn.xyz/ "Fengnan's Homepage") and [Tengyao
Wang](https://personal.lse.ac.uk/wangt60/ "Tengyao's Homepage"))

Implementation of the changepoint localization methods via the complementary 
sketching algorithm (collectively named 'charcoal') for high-dimensional 
regression coefficients, where the regression coefficients need not be 
individually sparse.

Including Algorithms 1, 2, 3 and 4 from Gao and Wang (2022).[^1]  
A function to generate linear regression samples with (multiple) changepoint
is also provided.

## Description of files

### R package

In `./R/` and `./man/` folders. Can be installed via `devtools::install_github('gaofengnan/charcoal')` in `R`.

## Reference

[^1]: Gao, F. and Wang, T. (2022) Sparse changepoint detection in 
high-dimensional linear regression. _arXiv preprint_, arXiv:2208.06326.
