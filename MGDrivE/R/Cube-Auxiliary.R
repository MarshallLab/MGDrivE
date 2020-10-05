###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   Auxiliary Functions for Cubes
#   Marshall Lab
#   March 2018
#
###############################################################################

#' MGDrivE: Inheritance Cube
#'
#' To model an arbitrary number of genotypes efficiently in the same mathematical
#' framework, we use a 3-dimensional array structure (cube) where each axis
#' represents the following information:
#'  * x: female adult mate genotype
#'	* y: male adult mate genotype
#'	* z: proportion of the offspring that inherits a given genotype (layer)
#'
#' The cube structure gives us the flexibility to apply tensor operations to the
#' elements within our equations, so that we can calculate the stratified population
#' dynamics rapidly; and within a readable, flexible computational framework.
#' This becomes apparent when we define the equation we use for the computation
#' of eggs laid at any given point in time:
#'
#'	\deqn{\overline{O(T_x)} = \sum_{j=1}^{n} \Bigg( \bigg( (\beta*\overline{s} * \overline{ \overline{Af_{[t-T_x]}}}) * \overline{\overline{\overline{Ih}}} \bigg) * \Lambda  \Bigg)^{\top}_{ij}}
#'
#' In this equation, the matrix containing the number of mated adult females
#' (\eqn{\overline{\overline{Af}}}) is multiplied element-wise with each one of
#' the layers containing the eggs genotypes proportions expected from this cross
#' (\eqn{\overline{\overline{\overline{Ih}}}}). The resulting matrix is then
#' multiplied by a binary 'viability mask' (\eqn{\Lambda}) that filters out
#' female-parent to offspring genetic combinations that are not viable due to
#' biological impediments (such as cytoplasmic incompatibility). The summation of
#' the transposed resulting matrix returns us the total fraction of eggs resulting
#' from all the male to female genotype crosses (\eqn{\overline{O(T_x)}}).
#'
#' Note: For inheritance operations to be consistent within the framework the
#' summation of each element in the z-axis (this is, the proportions of each one
#' of the offspring's genotypes) must be equal to one.
#'
#' @section Drive-specific Cubes:
#'
#' An inheritance cube in an array object that specifies inheritance probabilities
#' (offspring genotype probability) stratified by male and female parent genotypes.
#' MGDrivE provides the following cubes to model different gene drive systems:
#'  * \code{\link{cubeOneLocusTA}}: 1 Locus Maternal-Toxin/Zygotic-Antidote System
#'  * \code{\link{cubeTwoLocusTA}}: 2 Locus Maternal-Toxin/Zygotic-Antidote System
#'  * \code{\link{cubeClvR}}: 1 Locus Cleave and Rescue (ClvR)
#'  * \code{\link{cubeClvR2}}: 2 Locus Cleave and Rescue (ClvR)
#'  * \code{\link{cubeHoming1RA}}: Homing Drive with 1 Resistance Allele
#'  * \code{\link{cubeHomingDrive}}: CRISPR (Clustered Regularly Interspaced Short Palindromic Repeats) with 2 Resistance Allele
#'  * \code{\link{cubeECHACR}}: ERACR-CHACR
#'  * \code{\link{cubeECHACRX}}: ECHACR, X-Linked
#'  * \code{\link{cubeKillerRescue}}: Killer-Rescue System
#'  * \code{\link{cubeMEDEA}}: MEDEA (Maternal Effect Dominant Embryonic Arrest)
#'  * \code{\link{cubeReciprocalTranslocations}}: Reciprocal Translocation
#'  * \code{\link{cubeRIDL}}: RIDL (Release of Insects with Dominant Lethality)
#'  * \code{\link{cubeMendelian}}: Mendelian
#'  * \code{\link{cubeSplitDrive}}: Split CRISPR drive
#'  * \code{\link{cubeTGD}}: trans-complementing Gene Drive
#'  * \code{\link{cubeTGDX}}: trans-complementing Gene Drive, X-Linked
#'  * \code{\link{cubeWolbachia}}: Wolbachia
#'
#' @section Functions for Cubes:
#'
#' We provide one auxiliary function to operate on cube objects.
#'  * \code{\link{cube2csv}}: Export slices of a cube to .csv format
#'
#' @name MGDrivE-Cube
NULL
#> NULL

#' Export a Cube to .csv
#'
#' Export a cube as multiple .csv files (one for each genotype; slices of z-axis).
#' This function will create the directory if it doesn't exist. Files are stored
#' as slice_(z-slice)_(genotype).csv
#'
#' @param cube A cube object (see \code{\link{MGDrivE-Cube}} for options)
#' @param directory Directory to write .csv files to
#' @param digits Number of significant digits to retain in .csv output
#'
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#' # output directory
#' oPath <- "path/to/write/output"
#'
#' # setup inheritance cube for export, using Mendelian as the example
#' cube <- cubeMendelian()
#'
#' # write out
#' cube2csv(cube = cube, directory = oPath, digits = 3)
#' }
#'
#' @export
cube2csv <- function(cube, directory, digits = 3){

  # get variables for later
  dir = file.path(directory)
  dims = dim(cube$ih)
  names = dimnames(cube$ih)[[3]]

  # create directory if it doesn't exist
  if(!dir.exists(dir)){
    dir.create(dir)
  }

  # loop over slices
  for(i in 1:dims[3]){
    # file name for each slice
    fname = paste0("slice_",formatC(x=i,width=6,format="d",flag="0"),"_",names[i],".csv")

    # write output
    write.table(x = format(cube$ih[,,i], digits=digits, scientific = FALSE),
                file = file.path(dir,fname), append = FALSE, quote = FALSE,
                sep = ",", row.names = TRUE, col.names = NA)
  }
} # end function

#' Generate and Modify Default Genotype-specific Parameters
#'
#' This is an internal function for cubes.
#'
#' @param gtype character vector of genotypes
#' @param eta genotype-specific mating fitness, handles assortative mating as well
#' @param phi genotype-specific sex ratio at emergence
#' @param omega genotype-specific multiplicative modifier of adult mortality
#' @param xiF genotype-specific female pupatory success
#' @param xiM genotype-specific male pupatory success
#' @param s genotype-specific fractional reduction(increase) in fertility
#'
cubeModifiers <- function(gtype, eta = NULL, phi = NULL, omega = NULL,
                          xiF = NULL, xiM = NULL, s = NULL){

  # check all numeric arguments to have proper bounds
#  if(any(eta < 0)) stop("eta values must be positive [X>0]")
  if(any(phi > 1) || any(phi < 0)) stop("phi values must be between [0,1]")
  if(any(omega < 0)) stop("omega values must be positive [X>0]")
  if(any(xiF > 1) || any(xiF < 0)) stop("xiF values must be between [0,1]")
  if(any(xiM > 1) || any(xiM < 0)) stop("xiM values must be between [0,1]")
  if(any(s < 0)) stop("s values must be positive [X>0]")

  # add parameters in the right place
  set <- function(gtype,vector,vectorNew){
    if(is.null(vectorNew)){
      return(vector)
    } else {
      if(any(!names(vectorNew) %in% gtype)){
        stop("genotype(s) do not match genotypes in cube; please check names of input genotype-specific parameters")
      }
      if(length(vectorNew)==1L){
        if(is.null(names(vectorNew))){
          vector[1:length(vector)] = unname(vectorNew)
        } else {
          vector[names(vectorNew)] = vectorNew
        }
      } else {
        vector[names(vectorNew)] = vectorNew
      }
      return(vector)
    }
  }

  ## genotype-specific modifiers
  size = length(gtype)
#  etaN = setNames(object = rep.int(x = 1, times = size), nm = gtype)      # genotype-specific mating fitness
  phiN = setNames(object = rep.int(x = 0.5, times = size), nm = gtype)     # genotype-specific sex ratio at emergence
  omegaN = setNames(object = rep.int(x = 1, times = size), nm = gtype)    # genotype-specific multiplicative modifier of adult mortality
  xiFN = setNames(object = rep.int(x = 1, times = size), nm = gtype)      # genotype-specific female pupatory success
  xiMN = setNames(object = rep.int(x = 1, times = size), nm = gtype)      # genotype-specific male pupatory success
  sN = setNames(object = rep.int(x = 1, times = size), nm = gtype)        # genotype-specific fractional reduction(increase) in fertility

  # add the user parameters
#  etaO = set(gtype,etaN,eta)
  phiO = set(gtype,phiN,phi)
  omegaO = set(gtype,omegaN,omega)
  xiFO = set(gtype,xiFN,xiF)
  xiMO = set(gtype,xiMN,xiM)
  sO = set(gtype,sN,s)

  ##########
  # generate eta matrix
  ##########
  if(is.null(eta)){
    # default behaviour
    # This generates a mating matrix with all mating equal
    etaO <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype,gtype))

  } else if(all(lengths(eta) == 2) ) {
    # set mating weights over specific male genotypes and all female genotypes
    # takes input list of mate genotype, and weight

    # pull out values
    mGeno <- unlist(lapply(X = eta, FUN = "[[", 1))
    weight <- as.numeric(lapply(X = eta, FUN = "[[", 2))

    # check for correctly specified genotypes
    if(!all(mGeno %in% gtype)){
      stop(paste0("Genotypes specified for male mate is not
                  one of the genotypes in the specified inheritance pattern."))
    }
    if(any(weight < 0)){
      stop(paste0("Mating ability weights must be greater than or equal to zero."))
    }

    # setup/fill mating matrix
    etaO <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype,gtype))
    etaO[ ,mGeno] <- rep(x = weight, each = size)

  } else if(all(lengths(eta) == 3) ) {
    # mating ability, depends on genotype of male and female
    # takes input list of female genotype, mate genotype, and weight

    # pull out values
    fGeno <- unlist(lapply(X = eta, FUN = "[[", 1))
    mGeno <- unlist(lapply(X = eta, FUN = "[[", 2))
    weight <- as.numeric(lapply(X = eta, FUN = "[[", 3))

    # check for correctly specified genotypes
    if(!all(fGeno %in% gtype)){
      stop(paste0("Genotypes specified for female mate is not
                  one of the genotypes in the specified inheritance pattern."))
    }
    if(!all(mGeno %in% gtype)){
      stop(paste0("Genotypes specified for male mate is not
                  one of the genotypes in the specified inheritance pattern."))
    }
    if(any(weight < 0)){
      stop(paste0("Mating ability weights must be greater than or equal to zero."))
    }

    # setup/fill mating matrix
    etaO <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype,gtype))
    etaO[cbind(fGeno,mGeno)] <- weight

  } else {
    stop("eta is incorrectly formatted.\n
         Check the examples or set as NULL for default behaviour.")
  }

  # return named list
  return(
    list(
      eta=etaO,
      phi=phiO,
      omega=omegaO,
      xiF=xiFO,
      xiM=xiMO,
      s=sO
    )
  )
}
