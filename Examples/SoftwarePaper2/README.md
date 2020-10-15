# Simulation Files

This directory contains R scripts that run the simulations presented as examples in our main text. For more detailed and commented information regarding the use of the `MGDrivE2` package, please read the vignettes in the package.

  * Data: folder containing time-series data used for time-varying carrying capacities and adult mortality.
  * `sims.R`: a script that runs the simulation once and generates a plot, in R.
  * `sims-server.R`: a script which ran 100 stochastic simulations on our server and returns the sampled trajectories as .CSV files.