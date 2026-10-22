## Test environments

* local Ubuntu 24.04, R 4.3.3
* GitHub Actions: windows-latest, macos-latest, ubuntu-latest (devel, release, oldrel-1)
* win-builder (devel)   
* R-hub / mac builder (optional)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* Possibly misspelled words in DESCRIPTION: Antonucci, Crocetta, Irpino
  (authors' names), Wasserstein, barycenters, Cronbach's (technical terms).
  They are spelled correctly.

## Notes

* The package includes C++ code (via Rcpp); tests and examples run cleanly
  under valgrind (0 errors).
* The vignettes use the 'univer' dataset of the CRAN package 'CUB' (in
  Suggests); they are evaluated only if 'CUB' is available.
