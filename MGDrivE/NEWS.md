# MGDrivE 1.1.0

### Major Changes

* data.table moved to Imports, from Depends.
* Parameter names made consistent in several auxiliary functions:
  * `retrieveOutput()`
  * `calcQuantiles()`
  * `plotMGDrivESingle()`

### Minor Changes

* Several spelling errors addressed.
* Realistic landscape example changed to use the zero-inflated kernel, instead of the basic exponential kernel.
* Several plotting tweaks.
* Added links to the accompanying data analysis package, [MoNeT-MGDrivE](https://pypi.org/project/MoNeT-MGDrivE/), part of the [MoNeT](https://chipdelmal.github.io/MoNeT/) package.



# MGDrivE 1.5.0

### Major Changes

* Complete internal rebuild of MGDrivE implementation.
  * The underlying mathematics are the same, only the implementation has been changed.
  * Significant memory reductions.
  * Significant computational reduction.
  * All internal objects and functions have been updated.
  * Most important, the stochastic implementation had a bug in it, which has been resolved.
  * All releases are now numeric vectors/matrices indicating the genotype and release number; this is handled internally.
  * It is now possible to release mated females.
  * Male-mating ability is now female-genotype dependent, allowing for assortative mating.
  * `parameterizeMGDrivE()` takes several new parameters
    * `sampTime` indicates how often output is written from the simulation
    * `inheritanceCube` is now required to parameterize the initial genotype distributions
    * `LarPopRatio`, `AdPopRatio_F`, and `AdPopRatio_M` have been updated internally to reflect accurate default behavior and handle different user input to set them.
* data.table has been removed from the dependencies. This implies rebuilds of:
  * `splitOutput()`
  * `aggregateFemales()`
  * `calcQuantiles()`
* `splitOutput()` and `aggregateFemales()` have been fixed to properly use the `writeDir` parameter.

### Minor Changes

* All `verbose` options have been updated for consistency. The default is `TRUE`.
* All parameters that previously had to be vectors the same length as the number of patches have been updated to take a single number, implying that all parameters are the same for each patch, or as a vector, so each patch can be specified individually.
* Plotting functions have been updated to handle any sampling scheme (i.e., if output is not written every day). 
* `parameterizeMGDrivE()` had internal loops replaced with vectorized functions.
* Spelling errors and documentation inconsistencies were addressed.
* Citation was updated to reflect publication in [Methods in Ecology and Evolution](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13318).



# MGDrivE 1.6.0

### Major Changes

* Rebuild of the Migration Function
  * Migration is handled at the Network level now. This removes objects from the Patch class, making the package lighter and more efficient.
  * Migration is no longer Dirichlet distributed.
* New Inheritance Patterns
  * One and Two locus Cleave and Rescue (ClvR) constructs have been made available.
  * ERACR/eCHACR constructs have been made available.

### Minor Changes

* Spelling checked and errors corrected.
* Function links in the documentation have been updated.



# MGDrivE 1.6.1

### Major Changes

* Rebuild of `oneDay_mating_deterministic_Patch()` and `oneDay_mating_stochastic_Patch()`. They properly account for zero probabilities with males and the associated female deaths.
  * This fixes a bug introduced when `eta()` was turned into a matrix object, allowing selective mating.
  * `omega()` is supposed to be a matrix, and change for each node. That isn't implemented here, and I don't know why not.

### Minor Changes

* `calcOmega()` has been updated internally.
  * A closed-form expression was derived, and the `stats::uniroot()` dependency is removed.
  * This does not actually fix the issue found, which is in the bounds for `mu`.
  * Better bounds-checking was introduced, with a more clear error message, to warn users when they encounter this error.



# MGDrivE 1.6.2

### Major Changes

* None

### Minor Changes

* Fixed deposition issues in `cubeClvR2()` - Thank you Dr. Leftwich.
* Updated internal documentation - per CRAN requests.
* Exported several new inheritance patterns that have been in development. 
* Updated website documentation, switched to the current version of `pkgdown`.
* Updated vignette options to reduce package size. 
