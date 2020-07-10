# `robttest`: Robust t-test alternative for small samples

Standard inference about a scalar parameter estimated via GMM amounts to applying a t-test to a particular set of observations.
If the number of observations is not very large, then moderately heavy tails can lead to poor behaviour of the t-test.
This Stata package implements the method described in [Müller (2020)](http://www.princeton.edu/~umueller/heavymean.pdf), which combines extreme value theory for the smallest and largest observations with a normal approximation for the average of the remaining observations to construct a more robust alternative to the t-test.
This new test controls size more succesfully in small samples compared to existing methods.
Analytical results in the canonical inference for the mean problem demonstrate that the new test provides a refinement over the full sample t-test under more than two but less than three moments, while the bootstrapped t-test does not.


# References

Müller, Ulrich K (2020). "A More Robust t-Test". Working Paper. https://www.princeton.edu/~umueller/heavymean.pdf.
