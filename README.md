# `robttest`: Robust t-test for potentially heavy-tailed observations

Standard inference about a scalar parameter estimated via GMM amounts to applying a t-test to a particular set of observations.
If the number of observations is not very large, then moderately heavy tails can lead to poor behaviour of the t-test.
This Stata package implements the method described in [Müller (2020)](http://www.princeton.edu/~umueller/heavymean.pdf), which combines extreme value theory for the smallest and largest observations with a normal approximation for the average of the remaining observations to construct a more robust alternative to the t-test.
This new test controls size more succesfully in small samples compared to existing methods.
It is implemented as a postestimation command that can be used after the STATA commands "regress", "ivregress", "areg" , "logit" or "probit" as long as these are used with the standard error option "robust" or "cluster". 

# Install

`robttest` is not currently available from SSC. To install directly from this repository, you can copy and run the following lines in Stata:
```stata
// Remove program if it existed previously
cap ado uninstall robttest
// Install most up-to-date version
net install robttest, from("https://raw.githubusercontent.com/acarril/robttest/master/src")
```

# References

Müller, Ulrich K (2020). "A More Robust t-Test". Working Paper. https://www.princeton.edu/~umueller/heavymean.pdf.
