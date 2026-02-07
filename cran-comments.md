## R CMD check results
0 errors | 0 warnings | 0 note

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages


## Resubmission

This is a resubmission. In this version I have:

1. Fixed \dontrun{} usage:
   - Replaced \dontrun with \donttest for examples taking >5 seconds
   - Removed unnecessary \dontrun wrappers

2. Removed specific seed setting within functions:
   - Removed set.seed() calls from FIRM() and performance() functions
   - Added seed as optional parameter (default NULL) where needed

3. Added proper attribution for third-party code:
   - Added Gad Abraham as contributor (ctb) and copyright holder (cph) 
     in Authors@R field
   - Added Netherlands eScience Center as copyright holder (cph)
   - Created inst/COPYRIGHTS file with detailed modification notes
   - Updated Description field to acknowledge flashpca source

