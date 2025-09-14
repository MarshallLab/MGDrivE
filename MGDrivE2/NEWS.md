# MGDrivE2 2.1.0

New features:

1. Decoupled/Imperial sampling can now calculate malaria-induced severe disease and mortality. See `human_Imperial_ODE` function for implementation details.
2. `equilibrium_Imperial_decoupled` updated to calculate the total number of mosquitoes given human population parameters.

Bugfixes:

1. `human_Imperial_ODE` updated to correctly calculate the per-capita force of infection. Does not fail with suppression constructs anymore due to divide-by-zero error.
2. `base_female_Imperial` updated with correct mosquito density equations to begin the mosquito population at equilibrium, based on epi parameters.

# MGDrivE2 2.0.0

New features

1. Decoupled epi-mosquito dynamics now implemented with Imperial model of malaria transmission. See Vignettes for details.
2. `add_interventions` allows mosquito lifecycle parameters to vary in the presence of LLINs and IRS.

# MGDrivE2 1.2.0

New features:

1. Function `batch_migration` now takes an additional argument `stages` to allow simulation of batch migration where only one of egg, larval, or pupal stages experiences batch migration events.
2. Stuff about decoupled epi-mosquito dynamics

Bugfixes, miscellaneous:

1. If tracking events using `sim_trajectory_base_CSV`, write to file called "Tracking.csv" instead of "Events.csv" to avoid clashes with summary function `split_aggregate_CSV`.

# MGDrivE2 1.1.0

New features:

1. Update to include batch migration, see new file `R/auxiliary-batch.R` and user-facing function `batch_migration` which sets up the batch migration events. Functions `sim_trajectory_R` and `sim_trajectory_CSV` now take an additional argument `batch` which may be set to `NULL` if no batch migration is to be simulated.
2. Added function `get_shape` to get the shape parameter (number of bins) for Erlang developmental delays, given rate parameter `q` and coefficient of variation `cv`.

Bugfixes, miscellaneous:

1. In `sim_trajectory_base_CSV` and `sim_trajectory_base_R` the hazards are now correctly evaluated at the start of the time step rather than at the end; [times[i-1],times[i]). Both functions no longer need the (extraneous) argument `dt`.
2. Remove argument `t0` from `sim_trajectory_base_CSV` and `sim_trajectory_base_R`, all simulations are assumed to start at time 0. Change named argument `tt` to `tmax`.
3. In `sim_trajectory_base_CSV` and `sim_trajectory_base_R` testing for whether to preform a release is now based on a `<=` check, rather than equality.
4. Add argument `n` to function `spn_T_lifecycle_network` to allow simulation of multiple nodes with no migration between them. Useful for debugging as well as for simulation of landscapes where nodes are far enough apart that batch migration is the only way for mosquitoes to travel between nodes.

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
