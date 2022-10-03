# `scpc`: Spatial Correlation Robust Inference

This Stata package implements the C-SCPC method described in [Müller and Watson (2021)](http://www.princeton.edu/~umueller/SCPC.pdf) and SCPC method [Müller and Watson (2022)](http://www.princeton.edu/~umueller/SpatialRegression.pdf) for the construction of confidence intervals that account for many forms of spatial correlation.
It is implemented as a postestimation command that can be used after the STATA commands "regress", "ivregress", "areg" , "logit" or "probit" as long as these are used with the standard error option "robust" or "cluster". If "cluster" is chosen, then the method assumes that all observations in a cluster are at the same spatial location, and corrects for potential spatial correlation between clusters. 

# Install

`scpc` is not currently available from SSC. To install directly from this repository, you can copy and run the following lines in Stata:
```stata
// Remove program if it existed previously
cap ado uninstall scpc
// Install most up-to-date version
net install scpc, from("https://raw.githubusercontent.com/ukmueller/SCPC/master/src")
```

# References

Müller, Ulrich K and Mark W. Watson (2021). "Spatial Correlation Robust Inference". Working Paper. https://www.princeton.edu/~umueller/SHAR.pdf.

Müller, Ulrich K and Mark W. Watson (2022). "Spatial Correlation Robust Inference in Linear Regression and Panel Models". Working Paper. https://www.princeton.edu/~umueller/SptialRegression.pdf.
