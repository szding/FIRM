## R CMD check results
0 errors | 0 warnings | 0 note

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages


## Resubmission

This is a resubmission. In this version I have:

* Fixed \dontrun{} usage: replaced with \donttest for examples >5s
* Removed set.seed() calls from functions; added seed parameter
* Added proper attribution for flashpca third-party code:
  - Gad Abraham [ctb, cph] in Authors@R
  - Created inst/COPYRIGHTS with modification notes
* Updated DESCRIPTION: refined author roles, corrected GitHub URL

