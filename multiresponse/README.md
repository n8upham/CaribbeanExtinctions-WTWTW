Scripts for: Where the wild things were: intrinsic and extrinsic extinction predictors in the world’s most depleted mammal fauna

Samuel T. Turvey1, Clare Duncan1,2, Nathan S. Upham3,4,5, Xavier Harrison6, Liliana M. Dávalos7,8

R code to reproduce the multinomial results. Address code queries to liliana.davalos@stonybrook.edu

1) Both data files should be in the same folder as the scripts, or else modify to give file paths,

2) phylo_ultra_trees.r generates ultrametric trees that can generate the cov matrices for the phylogenetic structure of errors,

3) phylo_stan_multi_method.r runs 4 models with different multinomial links,

4) phylo_stan_multi_models.r runs cumulative link models with different sets of covariates and imputed vs. missing data, requires RData file from phylo_stan_multi_method.r,

5) phylo_stan_multi_models_aca.r runs adjacent category link models, requires RData file from phylo_stan_multi_method.r,

6) model_comparison.R takes all the resulting RData files, adds the loo and waic criterion to the models, and compares them when feasible.