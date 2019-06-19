########################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   MGDrivE Auxiliary
#   Marshall Lab
#   November 2017
#
########################################################################


########################################################################
# Calculate Parameter Values
########################################################################

#' Solve for Omega (additional genotype-specific mortality)
#'
#' Solves for root of equation of geometrically-distributed lifespan for value of omega.
#'
#' @param mu daily mortality probability (discrete-time hazard, called \code{muAd} in code)
#' @param lifespanReduction percent reduction in lifespan from average lifespan
#' (target average lifespan will be \eqn{\frac{1}{\mu_{Ad}} \times lifespanReduction})
#'
#' @importFrom stats uniroot
#'
#' @examples
#' # reduce lifespan by 10%
#' #  Example mu is an average for Aedes
#' newOmega <- getOmega(mu = 0.11, lifespanReduction = 0.90)
#'
#' @export
getOmega <- function(mu, lifespanReduction){

  if(lifespanReduction > 1 | lifespanReduction < 0){
    stop("parameter 'lifespanReduction' must be in interval [0,1]")
  }

  lifespan <- (1/mu) * lifespanReduction

  getOmega_f <- function(omega,mu,lifespan){
    (1/(1-omega+(omega*mu)))-lifespan
  }

  return(uniroot(f=getOmega_f,interval=c(0,1),lifespan=lifespan,mu=mu,maxiter=1e4)$root)
}

########################################################################
# Delete files in a directory
########################################################################

#' Erase all files in a directory
#'
#' Given a directory path, check it exists and if so delete all its contents.
#'
#' @param directory directory whose contents will be deleted
#' @param verbose Chatty? Default is TRUE
#'
#' @examples
#' # Path to directory, can tilde expand
#' myPath <- ""
#'
#' # Erase directory
#' #  No return value
#' eraseDirectory(directory = myPath)
#'
#' @export
eraseDirectory <- function(directory, verbose = TRUE){
  # check directory exists
  if(!dir.exists(directory)){
    if(verbose){cat("no such directory exists\n")}
    return(NULL)
  }
  dirFiles = list.files(path = directory)
  # begin deleting contents
  if(length(dirFiles)>0){
      for(i in dirFiles){
        if(verbose){cat("removing file: ",file.path(directory, i),"\n", sep = "")}
        #file.remove(file.path(directory, i))
        unlink(x = file.path(directory, i), recursive = TRUE)
      }
  }
  # end deleting contents
}

########################################################################
# Post-processing of Output
########################################################################

#' Split Output by Patch
#'
#' Split output into multiple files by patches.
#'
#' @param readDir Directory where output was written to
#' @param writeDir Directory to write output to. Default is readDir
#' @param remFile Remove original output? Default is TRUE
#' @param numCores How many cores to use when writing output. Default is 1
#' @param verbose Chatty? Default is TRUE
#'
#' @import data.table
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' \dontrun{
#' # This example assumes user has already run MGDrivE and generated output.
#' #  If that's untree, see vignette for complete example
#' fPath <- "path/to/data/containing/folder"
#' oPath <- "path/to/write/output"
#'
#' # split data by patch, keep original files
#' #  not return value
#' #  numCores uses data.table::setDTthreads internally
#' splitOutput(readDir = fPath, writeDir = oPath, remFile = FALSE, numCores = 1)
#'
#' # Alternatively, remove the original files and write new ones in their place
#' fPath <- "path/to/data/containing/folder"
#'
#' splitOutput(readDir = fPath, writeDir = NULL, remFile = TRUE, numCores = 1)
#' }
#'
#' @export
splitOutput <- function(readDir, writeDir=NULL, remFile=TRUE, numCores=1, verbose=TRUE){

  # Set data.table thread use.
  oldThread = data.table::getDTthreads(verbose = FALSE)
  data.table::setDTthreads(threads = numCores)

  # get all files in the read directory
  dirFiles = list.files(path = readDir, pattern = "^[MF]_.*\\.csv$", full.names = TRUE)

  # check write directory
  if(is.null(writeDir)){writeDir <- readDir}

  # initialize text, progress bar below
  if(verbose){cat("  Splitting", length(dirFiles), "files.\n")}
  if(remFile){
    if(verbose){cat("  Removing original files.\n")}
  } else {
    if(verbose){cat("  Not removing original files.\n\n")}
  }

  # loop over files
  for(selectFile in dirFiles){

    # Read in files
    fileIn = data.table::fread(input = selectFile, sep = ",",
                               header = TRUE, verbose = FALSE, showProgress = FALSE,
                               data.table = TRUE, logical01 = FALSE)

    # progress bar stuff
    if(verbose){
      pb = txtProgressBar(min = 0,max = length(unique(fileIn$Patch)),style = 3)
    }
    pbVal = 0

    # loop over patches, split and write out
    for(patch in unique(fileIn$Patch)){

      patchName = sub(pattern = ".csv",
                      replacement = paste0("_Patch",formatC(x = patch, width = 3, format = "d", flag = "0"),".csv"),
                      x = selectFile, fixed = TRUE)

      data.table::fwrite(x = fileIn[Patch==patch][ ,Patch:=NULL],
                         file = patchName, logical01 = FALSE,
                         showProgress = FALSE, verbose = FALSE)

      # some indication that it's working
      pbVal = pbVal+1
      if(verbose){setTxtProgressBar(pb = pb, value = pbVal)}
    }

    # check if removing original output
    if(remFile){file.remove(selectFile)}

  } # end loop over files

  # Reset data.table thread use.
  data.table::setDTthreads(threads = oldThread)
}

#' Aggregate Female Output by Genotype
#'
#' Aggregate over male mate genotype to convert female matrix output into vector output.
#'
#' @param readDir Directory to read input from
#' @param writeDir Directory to write output to. Default is readDir
#' @param genotypes Character vector of possible genotypes; found in \code{driveCube$genotypesID}
#' @param remFile Boolean flag to remove original (unaggregated) file
#' @param numCores Number of cores when reading/writing. Default is 1.
#' @param verbose Chatty? Default is TRUE
#'
#' @examples
#' \dontrun{
#' # This example assumes user has already run MGDrivE and generated output.
#' #  This also assumes that the user has already split output by patch
#' # See vignette for complete example
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
#' # numCores uses data.table::setDTthreads internally
#' aggregateFemales(readDir= fPath, writeDi = NULL, genotypes = cube$genotypesID,
#'                  remFile = FALSE, numCores = 1)
#' }
#'
#' @export
aggregateFemales <- function(readDir, writeDir=NULL, genotypes, remFile=FALSE,
                             numCores=1, verbose=TRUE){

  # Set data.table thread use.
  oldThread = data.table::getDTthreads(verbose = FALSE)
  data.table::setDTthreads(threads = numCores)

  #if the computer has the just-in-time compiler, it makes a huge difference in
  # string searches
  if(pcre_config()["JIT"]){options(PCRE_use_JIT = TRUE)}

  # check write directory
  if(is.null(writeDir)){writeDir <- readDir}

  #get female files
  femaleFiles = list.files(path = readDir, pattern = "^F_.*\\.csv$", full.names = TRUE)

  # initialize progress bar and text
  if(verbose){cat("  Aggregating", length(femaleFiles), "files.\n")}
  if(remFile){
    if(verbose){cat("  Removing original files.\n")}
  } else {
    if(verbose){cat("  Not removing original files.\n\n")}
  }

  if(verbose){
    pb = txtProgressBar(min = 0,max = length(femaleFiles),style = 3)
  }
  pbVal = 0

  #loop over files
  for(selectFile in femaleFiles){


    #read in file
    thisOutput =  data.table::fread(input = selectFile, sep = ",",
                                   header = TRUE, verbose = FALSE, showProgress = FALSE,
                                   data.table = TRUE, logical01 = FALSE)

    #get time/patch and remove from data.table
    timePatch = thisOutput[,"Time"] # this keeps object type too, useful later
    thisOutputNames = names(thisOutput)

    #aggregate output
    aggregateOut = vapply(X = genotypes,FUN = function(x){
      cols = grep(pattern = paste0("^", x), x = thisOutputNames, perl = TRUE)
      thisOutput[, rowSums(.SD), .SDcols = cols]
    },FUN.VALUE = numeric(nrow(thisOutput)))

    #name file and write out
    fileName = sub(pattern = "_",replacement = "_Aggregate_",x = selectFile, fixed = TRUE)
    data.table::fwrite(x = cbind(timePatch,aggregateOut),
                       file = fileName, logical01 = FALSE,
                       showProgress = FALSE, verbose = FALSE)

    # check if removing file
    if(remFile){file.remove(selectFile)}

    # some indication that it's working
    pbVal = pbVal+1
    if(verbose){setTxtProgressBar(pb = pb, value = pbVal)}

  }#end file loop

  # Reset data.table thread use.
  data.table::setDTthreads(threads = oldThread)

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
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' # This assumes user has run MGDrivE and output is in \code{readDir}
#' #  See vignette for examples on how to run MGDrivE
#'
#' # read/write dirs
#' fPath <- "folder/containing/output"
#' oPath <- "folder/to/write/stuff"
#'
#' # first, split output by patch and aggregate females by mate genotype
#' # remember, cube is for example and changes with simulation
#' cube <- cubeMendelian()
#'
#' splitOutput(readDir = fPath, writeDir = NULL, remFile = TRUE, numCores = 1)
#' aggregateFemales(readDir= fPath, writeDi = NULL, genotypes = cube$genotypesID,
#'                  remFile = FALSE, numCores = 1)
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


  # aggregate files
  for(sFile in 2:length(mFiles)){

    # males
    mMat <- mMat + matrix(data = scan(file = mFiles[sFile], what = integer(), sep = ",", skip = 1, quiet = TRUE),
                          ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))

    # females
    fMat <- fMat + matrix(data = scan(file = fFiles[sFile], what = integer(), sep = ",", skip = 1, quiet = TRUE),
                          ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))

  } # end loop over files

  # fix time addition
  mMat[ ,1] <- fMat[ ,1] <- 1:dim(mMat)[1]

  # print males
  fileName <- file.path(writeDir, file.path("M_LandscapeAgg_Run", hold[length(hold)], ".csv", fsep = ""))
  write.csv(x = mMat, file = fileName, row.names = FALSE)

  # print females
  fileName <- file.path(writeDir, file.path("F_LandscapeAgg_Run", hold[length(hold)], ".csv", fsep = ""))
  write.csv(x = fMat, file = fileName, row.names = FALSE)

} # end aggregate over landscape




#' Retrieve Output
#'
#' Read in output from directory. The resulting object will be a nested list;
#' outermost nesting dimension indexes runID, within runID elements are split by sex
#' and innermost nesting is over patches.
#'
#' @param readDir directory where output was written to; must not end in path seperator
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
#' Given the readDirectory, this function assumes the follow file structure: \cr
#'  * readDirectory
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
#' {M/F}_Mean_(patchNum).csv and {M/F}_Quantile_(quantNum)_(patchNum).csv.
#'
#' @param readDir Directory to find repetition folders in
#' @param writeDirectory Directory to write output
#' @param mean Boolean, calculate mean or not. Default is TRUE
#' @param quantiles Vector of quantiles to calculate. Default is NULL
#' @param numCores Number of cores when reading/writing. Default is 1.
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
#' calcQuantiles(readDir = fPath, writeDirectory = oPath, mean = TRUE,
#'               quantiles = NULL, numCores = 1)
#'
#' # here, calculate 2.5% and 97.5% quantiles
#' calcQuantiles(readDir = fPath, writeDirectory = oPath, mean = FALSE,
#'               quantiles = c(0.025, 0.975), numCores = 1)
#' }
#'
#' @return Writes output to files in writeDirectory
#' @export
calcQuantiles <- function(readDir, writeDirectory, mean=TRUE, quantiles=NULL,
                          numCores=1, verbose=TRUE){

  #safety check
  if(!mean && is.null(quantiles)){
    stop("User needs to specify the mean or which quantiles to calculate. ")
  }

  # Set data.table thread use.
  oldThread = data.table::getDTthreads(verbose = FALSE)
  data.table::setDTthreads(threads = numCores)

  #get files
  repFiles = list.dirs(path = readDir, full.names = TRUE, recursive = FALSE)
  patchFiles = lapply(X = repFiles, FUN = list.files, pattern = "^[FM].*\\.csv$")

  #subset females/males
  malePatches <- lapply(X = patchFiles, FUN = grep, pattern = "M_", fixed = TRUE, value=TRUE)
  femalePatches <- lapply(X = patchFiles, FUN = grep, pattern = "F_Aggregate", fixed = TRUE, value=TRUE)

  #generate a list of all patches to run over
  patchList = unique(regmatches(x = patchFiles[[1]],
                                m = regexpr(pattern = "Patch[0-9]+",
                                            text = patchFiles[[1]],
                                            perl = TRUE)))

  #read in a file initially to get variables and setup return array
  testFile <- data.table::fread(input = file.path(repFiles[1], patchFiles[[1]][1]),
                                header = TRUE, verbose = FALSE, showProgress = FALSE,
                                logical01 = FALSE, sep = ",", drop = "Time")

  #bunch of constants that get used several times
  numReps <- length(repFiles)
  columnNames <- c("Time", names(testFile))
  numRow <- dim(testFile)[1]
  numCol <- dim(testFile)[2]+1

  #setup input data holder
  popDataMale <- array(data = 0, dim = c(numRow,  numReps, numCol-1))
  popDataFemale <- array(data = 0, dim = c(numRow, numReps, numCol-1))

  #setup output data holder
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
  outputDataMale[,1,] <- 1:numRow
  outputDataFemale[,1,] <- 1:numRow


  # write initial outputs and setup progress bar
  if(verbose){cat("  Patches:", length(malePatches[[1]]), "\n")}
  if(verbose){cat("  Repetitions:", numReps, "\n\n")}

  if(verbose){
    pb = txtProgressBar(min = 0,max = length(malePatches),style = 3)
  }
  pbVal = 0


  #loop over all patches and do stats.
  for(patch in patchList){
    #get male and female files, all repetitions of this patch
    maleFiles <- vapply(X = malePatches,
                        FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                        FUN.VALUE = character(length = 1L))

    femaleFiles <- vapply(X = femalePatches,
                          FUN = grep, pattern = patch, fixed=TRUE, value=TRUE,
                          FUN.VALUE = character(length = 1L))


    #Read in all repetitions for this patch
    for(repetition in 1:numReps){
      popDataMale[ ,repetition, ] <- as.matrix(data.table::fread(input = file.path(repFiles[repetition], maleFiles[repetition]),
                                                                 header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                 logical01 = FALSE, sep = ",", drop = "Time"))
      popDataFemale[ ,repetition, ] <- as.matrix(data.table::fread(input = file.path(repFiles[repetition], femaleFiles[repetition]),
                                                                   header = TRUE, verbose = FALSE, showProgress = FALSE,
                                                                   logical01 = FALSE, sep = ",", drop = "Time"))
    }


    #if the user wants the mean, calculate and put out.
    if(mean){
      for(whichCol in 1:(numCol-1)){
        outputDataMale[ ,whichCol+1,1] <- .rowMeans(x = popDataMale[ , ,whichCol],
                                                           m = numRow, n = numReps)
        outputDataFemale[ ,whichCol+1,1] <- .rowMeans(x = popDataFemale[ , ,whichCol],
                                                             m = numRow, n = numReps)
      }

      #write output
      maleFileName <- file.path(writeDirectory,
                                file.path("M_Mean_", patch, ".csv", fsep = "")
      )
      femaleFileName <- file.path(writeDirectory,
                                  file.path("F_Mean_", patch, ".csv", fsep = "")
      )

      data.table::fwrite(x = as.data.frame(outputDataMale[ , ,1]),
                         file = maleFileName, verbose = FALSE,
                         showProgress = FALSE, logical01 = FALSE)
      data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,1]),
                         file = femaleFileName, verbose = FALSE,
                         showProgress = FALSE, logical01 = FALSE)
    }#end mean


    #if the user wants quantiles, do them and write output
    if(!is.null(quantiles)){
      for(whichCol in 1:(numCol-1)){
        outputDataMale[ ,whichCol+1, ] <- quantileC(Trials = popDataMale[ , ,whichCol],
                                                    Probs = quantiles)

        outputDataFemale[ ,whichCol+1, ] <- quantileC(Trials = popDataFemale[ , ,whichCol],
                                                      Probs = quantiles)

      }#end loop to calculate quantiles

      #write output
      for(whichQuant in 1:length(quantiles)){
        #file names
        maleFileName <- file.path(writeDirectory,
                                  file.path("M_Quantile_",
                                            formatC(x = quantiles[whichQuant], digits = 4,
                                                    format = "f", decimal.mark = "",
                                                    big.mark = NULL),
                                            "_", patch, ".csv", fsep = "")
        )
        femaleFileName <- file.path(writeDirectory,
                                    file.path("F_Quantile_",
                                              formatC(x = quantiles[whichQuant], digits = 4,
                                                      format = "f", decimal.mark = "",
                                                      big.mark = NULL),
                                              "_", patch, ".csv", fsep = "")
        )

        #write output
        data.table::fwrite(x = as.data.frame(outputDataMale[ , ,whichQuant]),
                           file = maleFileName, verbose = FALSE,
                           showProgress = FALSE, logical01 = FALSE)
        data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,whichQuant]),
                           file = femaleFileName, verbose = FALSE,
                           showProgress = FALSE, logical01 = FALSE)
      }#end loop to write files
    }#end quantiles

    # some indication that it's working
    pbVal = pbVal+1
    if(verbose){setTxtProgressBar(pb = pb, value = pbVal)}

  }#end loop over patches

  # Reset data.table thread use.
  data.table::setDTthreads(threads = oldThread)

}#end function

########################################################################
# Population Objects
########################################################################

#' Normalise a Numeric Vector
#'
#' Normalise a numeric vector to sum to one
#'
#' @param vector numeric vector
#'
normalise <- function(vector){
  if(all(vector==0)){
    return(vector)
  } else {
    return(vector/sum(vector))
  }
}

#' Create a Named Vector
#'
#' Create a named vector of 0s
#'
#' @param genotypesID character vector of possible genotypes
#'
createNamedPopVector <- function(genotypesID){
  out = rep(0,times=length(genotypesID))
  names(out) = genotypesID
  return(out)
}

#' Create a Named Matrix
#'
#' Create a named matrix of 0s
#'
#' @param genotypesID character vector of possible genotypes
#'
createNamedPopMatrix <- function(genotypesID){
  return(
    matrix(data=0,nrow=length(genotypesID),ncol=length(genotypesID),dimnames=list(genotypesID,genotypesID))
  )
}

#' Create a population array of vectors
#'
#' Creates an array for the population history to be stored. The length of the
#' array is equal to the window required for the model to run
#' (in our specific case it is equal to the sum of aquatic stages lengths).
#'
#' @param genotypesID character vector of possible genotypes
#' @param memoryWindow integer size of list structure
#'
initPopVectorArray <- function(genotypesID,memoryWindow){
  array=vector(mode="list",length=memoryWindow)
  for(i in 1:length(array)){array[[i]]=createNamedPopVector(genotypesID)}
  return(array)
}

#' Create a population array of matrices
#'
#' Creates an array for the population history to be stored. The length of the
#' array is equal to the window required for the model to run
#' (in our specific case it is equal to the sum of aquatic stages lengths).
#'
#' @param genotypesID character vector of possible genotypes
#' @param memoryWindow integer size of list structure
#'
initPopMatrixArray <- function(genotypesID,memoryWindow){
  array=vector(mode="list",length=memoryWindow)
  for(i in 1:length(array)){array[[i]]=createNamedPopMatrix(genotypesID)}
  return(array)
}

#' Create a primed population array of vectors
#'
#' Primes an array for the population history to be stored. The length of the
#' array is equal to the window required for the model to run
#' (in our specific case it is equal to the sum of aquatic stages lengths).
#'
#' @param primingVector a named vector population
#' @param memoryWindow integer size of list structure
#'
primePopVectorArray <- function(primingVector,memoryWindow){
  array=vector(mode="list",length=memoryWindow)
  for(i in 1:length(array)){array[[i]]=primingVector}
  return(array)
}

#' Create a primed population array of matrices
#'
#' Primes an array for the population history to be stored. The length of the
#' array is equal to the window required for the model to run
#' (in our specific case it is equal to the sum of aquatic stages lengths).
#'
#' @param primingMatrix a named matrix population
#' @param memoryWindow integer size of list structure
#'
primePopMatrixArray <- function(primingMatrix,memoryWindow){
  array=vector(mode="list",length=memoryWindow)
  for(i in 1:length(array)){array[[i]]=primingMatrix}
  return(array)
}

########################################################################
# Kernel-related
########################################################################

#' Calculates the zero-inflation part of a hurdle exponential kernel.
#'
#' Given the probability of an adult mosquito to stay in the same patch throughout its whole lifespan, and its mortality, it
#' calculates the height of the pulse-density part of the hurdle kernel.
#'
#' @param stayThroughLifespanProbability Probability of a mosquito to spend its whole lifespan in the same node
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
#' kernMat = calcHurdleExpKernel(distMat = distMat, rate = 10, pi = hHeight)
#'
#' @export
calcZeroInflation <- function(stayThroughLifespanProbability,adultMortality){
  stayThroughLifespanProbability^(adultMortality)
}
