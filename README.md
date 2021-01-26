# `scpc`: Spatial Correlation Robust Inference

This Stata package implements the SCPC method described in [Müller and Watson (2021)](http://www.princeton.edu/~umueller/SCPC.pdf) for the construction of confidence intervals that account for many forms of spatial correlation.
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

Müller, Ulrich K and Mark W. Watson (2021). "Spatial Correlation Robust Inference". Working Paper. https://www.princeton.edu/~umueller/SCPC.pdf.
