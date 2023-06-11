# Bayesian Record Linkage with Variables in One File

Accompanying illustrative code for the article "Bayesian Record Linkage with Variables in One File" by Gauri Kamat, Mingyang Shan, and Roee Gutman. 

Three files are provided corresponding to three general simulation scenarios: 
1. BRLVOF_correct.R: corresponds to the case when the true relationship between the covariate and the outcome is linear
2. BRLVOF_W.R: corresponds to the case when the true relationship between the covariate and the outcome includes an unobserved variable W
3. BRLVOF_nonlin.R: corresponds to the case when the true relationship between the covariate and the outcome is non-linear

The error level and regression parameters can be changed by resetting the parameters at the beginning of each R script. These parameters correspond to the values specified in Sizes.R.

Additional simulations can be constructed using the configurations outlined in the paper along with the code provided.
