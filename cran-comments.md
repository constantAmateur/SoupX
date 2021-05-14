## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Resubmission

This is a resubmission.  In this version I have:

* Corrected the issue that caused SoupX to be archived.  Specifically, I had not set the LazyDataCompression flag in DESCRIPTION to 'xz', causing the package to fail size checks.  This flag has been set.

* Removed examples from non-exported functions.

