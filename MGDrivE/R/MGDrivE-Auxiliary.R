###############################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   MGDrivE Auxiliary
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################
#######################################
# Calculate Parameter Values
#######################################

#' Solve for Omega (additional genotype-specific mortality)
#'
#' Solves for root of equation of geometrically-distributed lifespan for value of omega.
#'
#' @param mu Daily mortality probability (discrete-time hazard, called \code{muAd} in code)
#' @param lifespanReduction Target reduced lifespan, between 0 and 1
#' (target average lifespan will be \eqn{\frac{1}{\mu_{Ad}} \times lifespanReduction})
#'
#' @examples
#' # reduce lifespan by 10%
#' #  Example mu is an average for Aedes
#' newOmega <- calcOmega(mu = 0.11, lifespanReduction = 0.90)
#'
#' @export
calcOmega <- function(mu, lifespanReduction){

  ## Notes
  # The expected lifespan is 1/mu
  # The reduced lifespan is (1/mu)*lifespanReduction
  # Since mu is a daily death rate, and omega is a probability of life, to find
  #  the new omega, we need to use (1-mu)*omega.
  # Then, to convert back to our new expected lifespan, we have to subtract from
  #  1 again, so our new expected lifespan is 1/(1-(1-mu)*omega)
  # Thus, we solve (1/mu)*lifespanReduction = 1/(1-(1-mu)*omega) for omega

  # The mu lower bound is because of the series convergence used to derive
  #  1/mu as the expecte lifespan.
  # The upper bound may not exist.
  #  We require (1-mu)*omega <= 1, for probabilities to work within the simulation.
  #  Substituting our solution for omega and solving for lifespanReduction, we
  #  get that there is no bound.

  # Mapping lifespanReduction on [0,1) to [mu,inf) is an exercise left to the readers.

  # safety check
  if(lifespanReduction > 1 | lifespanReduction < mu){
    stop("parameter 'lifespanReduction' must be in interval [mu,1]")
  }

  # Calculate and return
  return( (lifespanReduction-mu)/(lifespanReduction*(1.0-mu)) )
}

#######################################
# Delete files in a directory
#######################################

#' Erase all files in a directory
#'
#' Given a directory path, check that it exists, and if so, delete all its contents.
#'
#' @param directory Directory whose contents will be deleted
#' @param verbose Chatty? Default is TRUE
#'
#' @examples
#' \dontrun{
#' # Path to directory, can tilde expand
#' myPath <- "~/path/to/write/output"
#'
#' # Erase directory
#' #  No return value
#' eraseDirectory(directory = myPath)
#' }
#'
#' @export
eraseDirectory <- function(directory, verbose = TRUE){

  # check directory exists
  if(!dir.exists(directory)){
    if(verbose){cat("no such directory exists\n")}
  }

  # list files/subdirectories to delete
  dirFiles = list.files(path = directory)

  # delete contents
  if(length(dirFiles)>0){
      for(i in dirFiles){
        if(verbose){cat("removing file: ",file.path(directory, i),"\n", sep = "")}
        unlink(x = file.path(directory, i), recursive = TRUE)
      }
  }

}

#######################################
# Population Objects
#######################################

#' Normalise a Numeric Vector
#'
#' Normalise a numeric vector to sum to one
#'
#' @param vector Numeric vector
#'
#' @keywords internal
#'
normalise <- function(vector){
  if(all(vector==0)){
    return(vector)
  } else {
    return(vector/sum(vector))
  }
}

#######################################
# Kernel-related
#######################################

#' Calculates the zero-inflation part of a hurdle exponential kernel.
#'
#' Given the probability of an adult mosquito to stay in the same patch throughout
#' its whole lifespan, and its mortality, it calculates the height of the pulse-density
#' part of the hurdle kernel.
#'
#' @param stayThroughLifespanProbability Probability of a mosquito to spend its
#' whole lifespan in the same node
#' @param adultMortality Adult mortality rate
#'
#' @examples
#' # setup distance matrix
#' # two-column matrix with latitude/longitude, in degrees
#' latLong = cbind(runif(n = 5, min = 0, max = 90),
#'                 runif(n = 5, min = 0, max = 180))
#'
#' # Vincenty Ellipsoid  distance formula
#' distMat = calcVinEll(latLongs = latLong)
#'
#' # get hurdle height
#' # Lets assume 80% stay probs and adult mortality of 0.1
#' hHeight <- calcZeroInflation(stayThroughLifespanProbability = 0.80,
#'                              adultMortality = 0.1)
#'
#' # calculate hurdle exponential distribution over distances
#' kernMat = calcHurdleExpKernel(distMat = distMat, rate = 10, p0 = hHeight)
#'
#' @export
calcZeroInflation <- function(stayThroughLifespanProbability,adultMortality){
  return(stayThroughLifespanProbability^(adultMortality))
}

#######################################
# Post-processing of Output
#######################################

#' Split Output by Patch
#'
#' Split output into multiple files by patches.
#'
#' @param readDir Directory where output was written to
#' @param writeDir Directory to write output to. Default is readDir
#' @param remFile Remove original output? Default is TRUE
#' @param verbose Chatty? Default is TRUE
#'
#' @importFrom utils setTxtProgressBar txtProgressBar write.table
#'
#' @examples
#' \dontrun{
#' # This example assumes user has already run MGDrivE and generated output.
#' #  If that's untree, see vignette for complete example
#' fPath <- "path/to/data/containing/folder"
#' oPath <- "path/to/write/output"
#'
#' # split data by patch, keep original files
#' #  no return value
#' splitOutput(readDir = fPath, writeDir = oPath, remFile = FALSE)
#'
#' # Alternatively, remove the original files and write new ones in their place
#' fPath <- "path/to/data/containing/folder"
#'
#' splitOutput(readDir = fPath, writeDir = NULL, remFile = TRUE)
#' }
#'
#' @export
splitOutput <- function(readDir, writeDir=NULL, remFile=TRUE, verbose=TRUE){

  # get all files in the read directory
  dirFiles = list.files(path = readDir, pattern = "^[MF]_.*\\.csv$", full.names = TRUE)

  # check write directory
  if(is.null(writeDir)){writeDir <- readDir}

  # initialize text, progress bar below
  if(verbose){
    # what we're doing
    cat("  Splitting", length(dirFiles), "files.\n")

    # whether or not we remove files
    if(remFile){
      cat("  Removing original files.\n")
    } else {
      cat("  Not removing original files.\n\n")
    }
  }

  # loop over files
  for(selectFile in dirFiles){

    # Read in files
    columnNames <- scan(file = selectFile, what = character(), sep = ",",
                        nlines = 1, quiet = TRUE)
    fileIn = matrix(data = scan(file = selectFile, what = numeric(), sep = ",",
                                skip = 1, quiet = TRUE),
                    ncol = length(columnNames), byrow = TRUE,
                    dimnames = list(NULL, columnNames))

    # get stuff for loop
    patches <- unique(fileIn[ ,2])
    patchRows <- seq.int(from = 1, to = dim(fileIn)[1], by = length(patches))
    fileBase <- strsplit(x = basename(path = selectFile), split = ".", fixed = TRUE)[[1]][1]

    # progress bar stuff
    if(verbose){
      pb = txtProgressBar(min = 0,max = length(patches),style = 3)
    }
    pbVal = 0

    # loop over patches, split and write out
    for(patch in patches){

      # output file name
      patchName = file.path(writeDir, paste0(fileBase, "_Patch",
                                             formatC(x = patch, width = 3, format = "d", flag = "0"),
                                             ".csv") )

      # write output
      #  subset matrix by patch number, then remove patch label column
      write.table(x = fileIn[patchRows, ][ ,-2], file = patchName, sep = ",",
                  row.names = FALSE, col.names = TRUE)

      # increment rows
      patchRows[] <- patchRows + 1L

      # some indication that it's working
      pbVal = pbVal+1
      if(verbose){setTxtProgressBar(pb = pb, value = pbVal)}
    }

    # check if removing original output
    if(remFile){file.remove(selectFile)}

  } # end loop over files

}

#' Aggregate Female Output by Genotype
#'
#' Aggregate over male mate genotype to convert female matrix output into vector output.
#'
#' @param readDir Directory to read input from
#' @param writeDir Directory to write output to. Default is readDir
#' @param genotypes Character vector of possible genotypes; found in \code{driveCube$genotypesID}
#' @param remFile Boolean flag to remove original (unaggregated) file
#' @param verbose Chatty? Default is TRUE
#'
#'
#' @examples
#' \dontrun{
#' # This example assumes user has already run MGDrivE and generated output.
#' #  This also assumes that the user has already split output by patch.
#' # See vignette for complete example.
#'
#' # set read/write directory
#' fPath <- "path/to/data/containing/folder"
#'
#' # Need genotypes from the cube run in the simulation
#' #  This is dependent on the simulation run
#' #  Using Mendelian cube for this example
#' cube <- cubeMendelian()
#'
#' # no return value from function
#' aggregateFemales(readDir= fPath, writeDir = NULL, genotypes = cube$genotypesID,
#'                  remFile = TRUE)
#' }
#'
#' @export
aggregateFemales <- function(readDir, writeDir=NULL, genotypes, remFile=TRUE,
                             verbose=TRUE){

  # check write directory
  if(is.null(writeDir)){writeDir <- readDir}

  #get female files
  femaleFiles = list.files(path = readDir, pattern = "^F_.*\\.csv$", full.names = TRUE)

  # initialize progress bar and text
  if(verbose){
    # what we're doing
    cat("  Aggregating", length(femaleFiles), "files.\n")

    # whether or not we remove shit
    if(remFile){
      cat("  Removing original files.\n")
    } else {
      cat("  Not removing original files.\n\n")
    }

    # progress bar
    pb = txtProgressBar(min = 0,max = length(femaleFiles),style = 3)
  }
  pbVal = 0


  # setup aggregation columns
  #  +1 offset for time column
  nGenos <- length(genotypes)
  idxList <- lapply(X = 0:(nGenos-1) * nGenos, FUN = "+", 2:(nGenos+1))


  #loop over files
  for(selectFile in femaleFiles){

    #read in file
    thisOutput = matrix(data = scan(file = selectFile, what = numeric(), sep = ",", skip = 1, quiet = TRUE),
                        ncol = 1+nGenos^2, byrow = TRUE)

    #aggregate output
    aggregateOut = vapply(X = idxList, FUN = function(x){
      rowSums(x = thisOutput[ ,x])
    }, FUN.VALUE = numeric(nrow(thisOutput)))

    #name file and write out
    fileName <- file.path(writeDir, sub(pattern = "_",replacement = "_Aggregate_",
                                        x = basename(path = selectFile), fixed = TRUE)
    )

    # combine time column from original file, with aggregated matrix of female counts
    write.table(x = cbind(thisOutput[ ,1],aggregateOut), file = fileName,
                sep = ",", row.names = FALSE, col.names = c("Time",genotypes))

    # check if removing file
    if(remFile){file.remove(selectFile)}

    # some indication that it's working
    pbVal = pbVal+1
    if(verbose){setTxtProgressBar(pb = pb, value = pbVal)}

  }#end file loop

}#end program


#' Aggregate Output Over Landscape
#'
#' This function aggregates the output of a run over the entire output, i.e., all
#' of the patches. It writes the output one level above the folder pointed to by
#' readDir, if writeDir is NULL. Output consists of 2 csv files, one for males and
#' one for females, "...M_LandscapeAgg_Run...csv".
#'
#' @usage aggregateOutput(readDir, writeDir=NULL)
#'
#' @param readDir Directory where output was written to
#' @param writeDir Directory to write output to. Default is one level above readDir
#'
#'
#' @examples
#' \dontrun{
#' # This assumes user has run MGDrivE and output is in fPath.
#' #  See vignette for examples on how to run MGDrivE
#'
#' # read/write dirs
#' fPath <- "folder/containing/output"
#' oPath <- "folder/to/write/stuff"
#'
#' # first, split output by patch and aggregate females by mate genotype
#' # remember, cube is for example and changes with simulation
#' #  landscape aggregation will work if females are not aggregated, but it's slower
#' cube <- cubeMendelian()
#'
#' splitOutput(readDir = fPath, writeDir = NULL, remFile = TRUE)
#' aggregateFemales(readDir= fPath, writeDi = NULL, genotypes = cube$genotypesID,
#'                  remFile = TRUE)
#'
#' # aggregate mosquitoes over entire landscape
#' #  no return value
#' aggregateOutput(readDir = fPath, writeDir = NULL)
#' }
#'
#' @export
aggregateOutput <- function(readDir, writeDir=NULL){

  # get all files in the read directory
  mFiles = list.files(path = readDir, pattern = "^M_.*\\.csv$", full.names = TRUE)
  fFiles = list.files(path = readDir, pattern = "^F_.*\\.csv$", full.names = TRUE)

  # check write directory
  hold <- strsplit(x = readDir, split = "/", fixed = TRUE)[[1]]
  if(is.null(writeDir)){
    writeDir <- paste0(hold[-length(hold)], collapse = "/")
  }

  # set object sizes for males and females
  columnNames <- scan(file = mFiles[1], what = character(), sep = ",",
                      quiet = TRUE, nlines = 1)

  mMat <- matrix(data = scan(file = mFiles[1], what = integer(), sep = ",", skip = 1, quiet = TRUE),
                 ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))


  columnNames <- scan(file = fFiles[1], what = character(), sep = ",",
                      quiet = TRUE, nlines = 1)

  fMat <- matrix(data = scan(file = fFiles[1], what = integer(), sep = ",", skip = 1, quiet = TRUE),
                 ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))

  # set things to reuse
  nRowM <- dim(mMat)[1]
  nColM <- dim(mMat)[2]
  nRowF <- dim(fMat)[1]
  nColF <- dim(fMat)[2]
  nReadM <- nRowM * nColM
  nReadF <- nRowF * nColF
  simTime <- mMat[ ,"Time"]

  # aggregate files
  for(sFile in 2:length(mFiles)){

    # males
    mMat <- mMat + matrix(data = scan(file = mFiles[sFile], what = integer(),
                                      sep = ",", skip = 1, n = nReadM, quiet = TRUE),
                          nrow = nRowM, ncol = nColM, byrow = TRUE)

    # females
    fMat <- fMat + matrix(data = scan(file = fFiles[sFile], what = integer(),
                                      sep = ",", skip = 1, n = nReadF, quiet = TRUE),
                          nrow = nRowF, ncol = nColF, byrow = TRUE)

  } # end loop over files

  # fix time addition
  mMat[ ,1] <- fMat[ ,1] <- simTime

  # print males
  fileName <- file.path(writeDir, file.path("M_LandscapeAgg_Run", hold[length(hold)], ".csv", fsep = ""))
  write.table(x = mMat, file = fileName, sep = ",", row.names = FALSE, col.names = TRUE)

  # print females
  fileName <- file.path(writeDir, file.path("F_LandscapeAgg_Run", hold[length(hold)], ".csv", fsep = ""))
  write.table(x = fMat, file = fileName, sep = ",", row.names = FALSE, col.names = TRUE)

} # end aggregate over landscape


#' Retrieve Output
#'
#' Read in output from directory. The resulting object will be a nested list;
#' outermost nesting dimension indexes runID, within runID elements are split by sex
#' and innermost nesting is over patches.
#'
#' @param readDir Directory where output was written to; must not end in path separator
#' @param verbose Chatty? Default is TRUE
#'
#' @importFrom utils read.csv
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' # Example assumes user has run and analyzed MGDrivE.
#' #  See vignette for examples of how to do that.
#'
#' # set read directory
#' fPath <- "path/to/split/aggregated/output"
#'
#' # read in data as nested lists
#' dataList <- retrieveOutput(readDir = fPath)
#' }
#'
#' @return Nested List
#'
#' @export
retrieveOutput <- function(readDir, verbose=TRUE){
  dirFiles = list.files(path = readDir)

  runID = strsplit(x = dirFiles,split = "_")
  runID = unique(x = grep( pattern = "Run", x = unlist(x = runID), value = TRUE))

  output = setNames(object = vector(mode = "list",length = length(runID)), runID)
  for(run in runID){
    if(verbose){cat("processing ",run,"\n",sep="")}

    runFiles = dirFiles[grep(pattern = run,x = dirFiles)]
    patches = regmatches(x = runFiles, m = regexpr(pattern = "Patch[0-9]+", text = runFiles))
    patches = sort(unique(patches))

    output[[run]]$F = setNames(object = vector(mode = "list",length = length(patches)),
                               nm = patches)
    output[[run]]$M = setNames(object = vector(mode = "list",length = length(patches)),
                               nm = patches)

    # retrieve males for this run
    males = runFiles[grep(pattern = "M_",x = runFiles, fixed = TRUE, useBytes = TRUE)]
    for(male in males){
      thisPatch = regmatches(x = male, m = regexpr(pattern = "Patch[0-9]+", text = male))
      thisOutput = read.csv(file = file.path(readDir, male))
      thisOutput = as.matrix(thisOutput)
      rownames(thisOutput) = NULL
      output[[run]]$M[[thisPatch]] = thisOutput[,-1]
    }

    # retrieve females for this run
    females = runFiles[grep(pattern = "F_Aggregate",x = runFiles)]
    for(female in females){
      thisPatch = regmatches(x = female, m = regexpr(pattern = "Patch[0-9]+", text = female))
      thisOutput = read.csv(file = file.path(readDir, female))
      thisOutput = thisOutput[,-1]
      output[[run]]$F[[thisPatch]] = as.matrix(thisOutput)
    }

  }
  return(output)
}

#' Summary Statistics for Stochastic MGDrivE
#'
#' This function reads in all repetitions for each patch and calculates either
#' the mean, quantiles, or both. User chooses the quantiles, up to 4 decimal places,
#' and enters them as a vector. Quantiles are calculated empirically. (order does not matter)  \cr
#'
#' Given the readDir, this function assumes the follow file structure: \cr
#'  * readDir
#'    * repetition 1
#'      * patch 1
#'      * patch 2
#'      * patch 3
#'    * repetition 2
#'      * patch 1
#'      * patch 2
#'      * patch 3
#'    * repetition 3
#'    * repetition 4
#'    * ... \cr
#'
#' Output files are *.csv contain the mean or quantile in the file name, i.e.
#' \{M/F\}_Mean_(patchNum).csv and \{M/F\}_Quantile_(quantNum)_(patchNum).csv.
#'
#' @param readDir Directory to find repetition folders in
#' @param writeDir Directory to write output
#' @param mean Boolean, calculate mean or not. Default is TRUE
#' @param quantiles Vector of quantiles to calculate. Default is NULL
#' @param verbose Chatty? Default is TRUE
#'
#' @examples
#' \dontrun{
#' # This function assumes network$multRun() has been performed, or several
#' #  network$oneRun() have been performed and all of the data has been split
#' #  and aggregated.
#'
#' # read/write paths
#' fPath <- "path/to/folder/ofFolders/with/data"
#' oPath <- "my/path/output"
#'
#' # here, only calculate mean, no quantiles
#' #  no return value
#' calcQuantiles(readDir = fPath, writeDir = oPath, mean = TRUE,
#'               quantiles = NULL)
#'
#' # here, calculate 2.5% and 97.5% quantiles
#' calcQuantiles(readDir = fPath, writeDir = oPath, mean = FALSE,
#'               quantiles = c(0.025, 0.975))
#' }
#'
#' @return Writes output to files in writeDir
#' @export
calcQuantiles <- function(readDir, writeDir, mean=TRUE, quantiles=NULL,
                          verbose=TRUE){

  # safety check
  if(!mean && is.null(quantiles)){
    stop("User needs to specify the mean or which quantiles to calculate. ")
  }

  # get all files
  patchFiles <- lapply(X = list.dirs(path = readDir, full.names = TRUE, recursive = FALSE),
                       FUN = list.files, pattern = "^[FM].*\\.csv$", full.names = TRUE)

  #s ubset females/males
  malePatches <- lapply(X = patchFiles, FUN = grep, pattern = "M_", fixed = TRUE, value=TRUE)
  femalePatches <- lapply(X = patchFiles, FUN = grep, pattern = "F_Aggregate", fixed = TRUE, value=TRUE)

  # generate a list of all patches to run over
  patchList = unique(regmatches(x = patchFiles[[1]],
                                m = regexpr(pattern = "Patch[0-9]+",
                                            text = patchFiles[[1]],
                                            perl = TRUE)))

  # read in a file initially to get variables and setup return array
  testFile <- read.csv(file = patchFiles[[1]][1], header = TRUE)

  # bunch of constants that get used several times
  numReps <- length(patchFiles)
  columnNames <- names(testFile)
  numRow <- dim(testFile)[1]
  numCol <- dim(testFile)[2]
  numRead <- numRow*numCol

  # setup input data holder
  popDataMale <- array(data = 0, dim = c(numRow,  numReps, numCol-1))
  popDataFemale <- array(data = 0, dim = c(numRow, numReps, numCol-1))

  # setup output data holder
  if(is.null(quantiles)){
    outputDataMale <- array(data = 0, dim = c(numRow, numCol, 1),
                            dimnames = list(NULL, columnNames, NULL))
    outputDataFemale <- array(data = 0, dim = c(numRow, numCol, 1),
                              dimnames = list(NULL, columnNames, NULL))
  } else {
    outputDataMale <- array(data = 0, dim = c(numRow, numCol, length(quantiles)),
                            dimnames = list(NULL, columnNames, NULL))
    outputDataFemale <- array(data = 0, dim = c(numRow, numCol, length(quantiles)),
                              dimnames = list(NULL, columnNames, NULL))
  }

  # set time from testFile
  outputDataMale[,1,] <- testFile[ ,"Time"]
  outputDataFemale[,1,] <- testFile[ ,"Time"]

  # write initial outputs and setup progress bar
  if(verbose){
    # how many patches and repetitions are we running through
    cat("  Patches:", length(patchList), "\n  Repetitions:", numReps, "\n\n")

    # setup progress bar
    pb = txtProgressBar(min = 0,max = length(patchList),style = 3)
  }
  pbVal = 0


  # loop over all patches and do stats.
  for(patch in patchList){
    # get male and female files, all repetitions of this patch
    maleFiles <- vapply(X = malePatches,
                        FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                        FUN.VALUE = character(length = 1L))

    femaleFiles <- vapply(X = femalePatches,
                          FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                          FUN.VALUE = character(length = 1L))


    # read in all repetitions for this patch
    #  remove time column after reading in.
    for(repetition in 1:numReps){
      popDataMale[ ,repetition, ] <- matrix(data = scan(file = maleFiles[repetition],
                                                        what = numeric(), n = numRead,
                                                        sep = ",", skip = 1, quiet = TRUE),
                                            nrow = numRow, ncol = numCol, byrow = TRUE)[ ,-1]

      popDataFemale[ ,repetition, ] <- matrix(data = scan(file = femaleFiles[repetition],
                                                          what = numeric(), n = numRead,
                                                          sep = ",", skip = 1, quiet = TRUE),
                                              nrow = numRow, ncol = numCol, byrow = TRUE)[ ,-1]
    }

    # if the user wants the mean, calculate and put out.
    if(mean){
      for(whichCol in 1:(numCol-1)){
        outputDataMale[ ,whichCol+1,1] <- .rowMeans(x = popDataMale[ , ,whichCol],
                                                    m = numRow, n = numReps)
        outputDataFemale[ ,whichCol+1,1] <- .rowMeans(x = popDataFemale[ , ,whichCol],
                                                      m = numRow, n = numReps)
      }

      # write output
      maleFileName <- file.path(writeDir,
                                file.path("M_Mean_", patch, ".csv", fsep = "") )
      femaleFileName <- file.path(writeDir,
                                  file.path("F_Mean_", patch, ".csv", fsep = "") )

      write.table(x = outputDataMale[ , ,1], file = maleFileName, sep = ",",
                  row.names = FALSE, col.names =TRUE)
      write.table(x = outputDataFemale[ , ,1], file = femaleFileName, sep = ",",
                  row.names = FALSE, col.names =TRUE)
    } # end mean


    # if the user wants quantiles, do them and write output
    if(!is.null(quantiles)){
      for(whichCol in 1:(numCol-1)){
        outputDataMale[ ,whichCol+1, ] <- quantileC(Trials = popDataMale[ , ,whichCol],
                                                    Probs = quantiles)

        outputDataFemale[ ,whichCol+1, ] <- quantileC(Trials = popDataFemale[ , ,whichCol],
                                                      Probs = quantiles)

      }#end loop to calculate quantiles

      # write output
      for(whichQuant in 1:length(quantiles)){
        # file names
        maleFileName <- file.path(writeDir,
                                  file.path("M_Quantile_",
                                            formatC(x = quantiles[whichQuant], digits = 4,
                                                    format = "f", decimal.mark = "",
                                                    big.mark = NULL),
                                            "_", patch, ".csv", fsep = "") )
        femaleFileName <- file.path(writeDir,
                                    file.path("F_Quantile_",
                                              formatC(x = quantiles[whichQuant], digits = 4,
                                                      format = "f", decimal.mark = "",
                                                      big.mark = NULL),
                                              "_", patch, ".csv", fsep = "") )

        # write output
        write.table(x = outputDataMale[ , ,whichQuant], file = maleFileName, sep = ",",
                    row.names = FALSE, col.names =TRUE)
        write.table(x = outputDataFemale[ , ,whichQuant], file = femaleFileName, sep = ",",
                    row.names = FALSE, col.names =TRUE)
      } # end loop to write files
    } # end quantiles

    # some indication that it's working
    pbVal = pbVal+1
    if(verbose){setTxtProgressBar(pb = pb, value = pbVal)}

  } # end loop over patches

} # end function
