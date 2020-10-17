# MGDrivE2 1.0.1

Small fixes to address CRAN comments:

  1. `R/hazard-functions-mosy.R`: function `make_move_male_haz` has a fix such that the origin/destination
  probability is correctly calculated.
  2. `R/equilibrium-lifecycle.R`: function `basic_eq_life` modified to not use the `<<-` global assignment operator.
  3. `vignettes/output-storage.Rmd:` replaced absolute path for files with `tempdir()` to conform to
  the CRAN standard.
  4. Edit `DESCRIPTION` to conform to CRAN standard.
  5. Remove `dontrun` code from examples where not necessary and direct users to appropriate vignettes,
  otherwise run examples.

# MGDrivE2 1.0.0

Initial release & CRAN submission.
