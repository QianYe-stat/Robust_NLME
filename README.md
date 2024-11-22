# Robust_NLME

## This repo contains R code that was used for the real data analysis and simulations in my paper ``Jointly Modeling Means and Variances for Nonlinear Mixed Effects Models with Measurement Errors and Outliers''

-----------------------------------------------------

### What's in the repo

+ `/src`: source code 
  + `00_dependencies.R`: the dependencies packages  
  + `Rnlme.R`: the function which jointly models mean and variance for NLME. Measurement error model can be added if needed.
  + `get_sd_bootstrap1.R`: the function produces adjusted SD by using parametric bootstrapping, when the random effect $a_i$ (see Model (2) in the paper) follows a normal distribution
  + `get_sd_bootstrap2.R`: the function produces adjusted SD by using parametric bootstrapping, when the random effect $a_i$ (see Model (2) in the paper) follows an inverse $\chi^2$ distribution
  + `get_Jloglike.R`, `est_fixed.R`, etc: all the other R files are self-defined functions sourced in the above functions and main analysis. 
  
+ `/simulation`: R files for different simulation settings (*see details about the simulation design in my paper Section 5*)
  + `Setting_I`: Setting I
  + `Setting_II_df=3`: Setting II with df=3
  + `Setting_II_df=5`: Setting II with df=5
  + `Setting_III`: Setting III (*The corresponding results shown in the Supplementery Materials*)
  + `Setting_IV`: Setting IV (*The corresponding results shown in the Supplementery Materials*)
  
+ `/real data analysis`: R files for real data analysis (*see details in my paper Section 4*)
  +
  +
  +
  
+ `/data`: 
  + `/toy data`:
  + `/DTP data`: (**EMPTY**) DTP data used in the real analysis. The data that support the findings of this study are available on request from the corresponding author. The data are not publicly available due to privacy or ethical restrictions.
  
  
  

















