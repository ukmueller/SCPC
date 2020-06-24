Meeting with Ulrich
June 24, 2020

Main goals:
1. Produce a proper package: help file, error checks, etc.
    - MATA section probably doesn't need a lot of work, but feel free to modify Stata wrapper.
2. Some version of reghdfe 
    - Issue: `reghdfe` doesn't produce score. However, covariance matrix must be somewhere (outer product)

Misc:
- There is some confusion regarding namespaces.
    - Do they affect other programs? (probably not)
    - What about between Mata and Stata? (not sure)
- Are we worried about errors in the code?
    - Probably able to upload a preliminary version to Github before SSC
    - Even in SSC you can quickly fix bugs
- There is some ancillary file (forgot the name), we need to include that in the package
