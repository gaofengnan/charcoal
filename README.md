# charcoal - Novel changepoint detection algorithms in high-dimensional linear regression

(By [Fengnan Gao](https://gaofn.xyz/ "Fengnan's Homepage") and [Tengyao
Wang](https://personal.lse.ac.uk/wangt60/ "Tengyao's Homepage"))

Implementation of the changepoint localization methods via the complementary sketching algorithm (collectively named '**charcoal**') for high-dimensional regression coefficients, where the regression coefficients need not be individually sparse.

* Including Algorithms 1, 2, 3 and 4 from Gao and Wang (2022).[^1]  
* A function to generate linear regression samples with (multiple) changepoint(s) is also provided.

## Description of files

* The folder `singlecell` contains the real data example in the paper, which came from [Suo et al, 2022](https://doi.org/10.1038/s41587-023-01734-7 "Suo et al, 2022) and was distributed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/version4/).

### R package

In `./R/` and `./man/` folders. Can be installed via `devtools::install_github('gaofengnan/charcoal')` in `R`.

[^1]: Gao, F. and Wang, T. (2022) Sparse change detection in 
high-dimensional linear regression. *arXiv preprint*, [arXiv:2208.06326](https://arxiv.org/abs/2208.06326).
