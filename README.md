# Robust_NLME

## This repo contains R code that was used for the real data analysis and simulations in my paper ``Jointly Modeling Means and Variances for Nonlinear Mixed Effects Models with Measurement Errors and Outliers''

-----------------------------------------------------

### What's in the repo

+ `/src`: Source code 
  + `00_dependencies.R`: the dependencies packages  
  + `Rnlme.R`: the function which jointly models mean and variance for NLME. Measurement error model can be added if needed.
  + `get_sd_bootstrap1.R`: the function produces adjusted SD by using parametric bootstrapping, when the random effect $a_i$ (see Model (2) in the paper) follows a normal distribution
  + `get_sd_bootstrap2.R`: the function produces adjusted SD by using parametric bootstrapping, when the random effect $a_i$ (see Model (2) in the paper) follows an inverse $\chi^2$ distribution
  + `get_Jloglike.R`, `est_fixed.R`, etc: all the other R files are self-defined functions sourced in the above functions and main analysis. 
  
+ `/simulation`: files for different simulation settings (*see details about the simulation design in my paper*)
  + `Setting_I`: simulation for the Setting I
  + `Setting_II`: simulation for Setting II
  

















