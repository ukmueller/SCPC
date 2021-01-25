# `SHAR`: Spatial Correlation Robust Inference

This Stata package implements the method described in [Müller and Watson (2021)](http://www.princeton.edu/~umueller/SHAR.pdf) for the construction of confidence intervals that account for many forms of spatial correlation.
It is implemented as a postestimation command that can be used after the STATA commands "regress", "ivregress", "areg" , "logit" or "probit" as long as these are used with the standard error option "robust" or "cluster". If "cluster" is chosen, then the method assumes that all observations in a cluster is at the same spatial location, and corrects for the spatial correlation between clusters. 

# Install

`SHAR` is not currently available from SSC. To install directly from this repository, you can copy and run the following lines in Stata:
```stata
// Remove program if it existed previously
cap ado uninstall SHAR
// Install most up-to-date version
net install SHAR, from("https://raw.githubusercontent.com/ukmueller/SHAR/master/src")
```

# References

Müller, Ulrich K and Mark W. Watson (2021). "Spatial Correlation Robust Inference". Working Paper. https://www.princeton.edu/~umueller/SHAR.pdf.
