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
# Delete files in a directory
########################################################################

#' Erase all files in a directory
#'
#' Given a directory path, check it exists and if so delete all its contents.
#'
#' @param directory directory whose contents will be deleted
#'
#' @export
eraseDirectory <- function(directory){
  # check directory exists
  if(!dir.exists(directory)){
    cat("no such directory exists\n")
    return(NULL)
  }
  dirFiles = list.files(path = directory)
  # begin deleting contents
  if(length(dirFiles)>0){
      for(i in dirFiles){
        cat("removing file: ",file.path(directory, i),"\n", sep = "")
        file.remove(file.path(directory, i))
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
#' @param directory Directory where output was written to; must not end in path seperator
#' @param multiCore Write output using multiple cores? Default is FALSE
#'
#' @export
splitOutput <- function(directory, multiCore=FALSE){
  dirFiles = list.files(path = directory, pattern = ".*\\.csv$")

  if(multiCore){
    nThread <- getDTthreads()
  } else{
    nThread <- 1
  }

  for(selectFile in dirFiles){
    cat("processing ",selectFile,"\n",sep="")
    fileIn = data.table::fread(input = file.path(directory, selectFile))
    # for each file, get all the patches and split into multiple files

    for(patch in unique(fileIn$Patch)){
      patchIn = fileIn[Patch==patch]
      patchName = sub(pattern = ".csv",
                      replacement = paste0("_Patch",formatC(x = patch, width = 4, format = "d", flag = "0"),".csv"),
                      x = selectFile, fixed = TRUE)
      data.table::fwrite(x = patchIn, file = file.path(directory, patchName), nThread = nThread)
    }
    cat("removing ",selectFile,"\n",sep="")
    file.remove(file.path(directory, selectFile))
  }
}

#' Aggregate Female Output by Genotype
#'
#' Aggregate over male mate genotype to convert female matrix output into vector output.
#'
#' @param directory Directory where output was written to; must not end in path seperator
#' @param genotypes Character vector of possible genotypes; found in \code{driveCube$genotypesID}
#' @param remove Boolean flag to remove original (unaggregated) file
#' @param multiCore Write output using multiple cores? Default is FALSE
#'
#' @export
aggregateFemales <- function(directory, genotypes, remove=FALSE, multiCore=FALSE){

  #if the computer has the just-in-time compiler, it makes a huge difference in
  # string searches
  if(pcre_config()["JIT"]){options(PCRE_use_JIT = TRUE)}

  if(multiCore){
    nThread <- getDTthreads()
  } else{
    nThread <- 1
  }

  #get files and subset females
  dirFiles = list.files(path = directory, pattern = ".*\\.csv$")
  femaleFiles = dirFiles[grep(pattern = "AF1",fixed = TRUE,x = dirFiles)]

  #loop over files
  for(selectFile in femaleFiles){

    #read in file
    cat("processing ",selectFile,"\n",sep="")
    thisPatch = regmatches(x = selectFile, m = regexpr(pattern = "Patch[0-9]+", text = selectFile, perl = TRUE))
    thisOutput = data.table::fread(file = file.path(directory, selectFile))

    #get time/patch and remove from data.table
    timePatch = thisOutput[,c("Time", "Patch")]
    thisOutput = thisOutput[,-c("Time", "Patch")]
    thisOutputNames = names(thisOutput)

    #aggregate output
    aggregateOut = vapply(X = genotypes,FUN = function(x){
      cols = grep(pattern = paste0("^", x), x = thisOutputNames, perl = TRUE)
      thisOutput[, rowSums(.SD), .SDcols = cols]
    },FUN.VALUE = numeric(nrow(thisOutput)))

    #name file and write out
    fileName = sub(pattern = "_",replacement = "_Aggregate_",x = selectFile, fixed = TRUE)
    cat("writing ",fileName,"\n",sep="")
    data.table::fwrite(x = cbind(timePatch,aggregateOut),
                       file = file.path(directory, fileName), nThread = nThread)

    if(remove){
      cat("removing ",selectFile,"\n",sep="")
      file.remove(file.path(directory, selectFile))
    }#end if
  }#end file loop
}#end program

#' Retrieve Output
#'
#' Read in output from directory. The resulting object will be a nested list;
#' outermost nesting dimension indexes runID, within runID elements are split by sex
#' and innermost nesting is over patches.
#'
#' @param directory directory where output was written to; must not end in path seperator
#' @param genotypes character vector of possible genotypes; found in \code{driveCube$genotypesID}
#'
#' @export
retrieveOutput <- function(directory, genotypes){
  dirFiles = list.files(path = directory)

  runID = strsplit(x = dirFiles,split = "_")
  runID = unique(x = grep( pattern = "Run", x = unlist(x = runID), value = TRUE))

  output = setNames(object = vector(mode = "list",length = length(runID)), runID)
  for(run in runID){
    cat("processing ",run,"\n",sep="")

    runFiles = dirFiles[grep(pattern = run,x = dirFiles)]
    patches = regmatches(x = runFiles, m = regexpr(pattern = "Patch[0-9]+", text = runFiles))
    patches = sort(unique(patches))

    output[[run]]$F = setNames(object = vector(mode = "list",length = length(patches)),
                               nm = patches)
    output[[run]]$M = setNames(object = vector(mode = "list",length = length(patches)),
                               nm = patches)

    # retrieve males for this run
    males = runFiles[grep(pattern = "ADM",x = runFiles)]
    for(male in males){
      thisPatch = regmatches(x = male, m = regexpr(pattern = "Patch[0-9]+", text = male))
      thisOutput = read.csv(file = file.path(directory, male))
      thisOutput = as.matrix(thisOutput)
      rownames(thisOutput) = NULL
      output[[run]]$M[[thisPatch]] = thisOutput[,-c(1,2)]
    }

    # retrieve females for this run
    females = runFiles[grep(pattern = "AF1_Aggregate",x = runFiles)]
    for(female in females){
      thisPatch = regmatches(x = female, m = regexpr(pattern = "Patch[0-9]+", text = female))
      thisOutput = read.csv(file = file.path(directory, female))
      thisOutput = thisOutput[,-c(1,2)]
      output[[run]]$F[[thisPatch]] = as.matrix(thisOutput)
    }

  }
  return(output)
}

#' Summary Statistics for Stochastic MGDrivE
#'
#' This function reads in all repetitions for each patch and calculates either
#' the mean, quantiles, or both. User chooses the quantiles, up to 4 decimal places,
#' and enters them as a vector. (order does not matter)  \cr
#'
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
#'    * ...
#'
#' @param readDirectory Directory to find repetition folders in
#' @param writeDirectory Directory to write output
#' @param mean Boolean, calculate mean or not. Default is TRUE
#' @param quantiles Vector of quantiles to calculate. Default is NULL
#'
#' @return Writes output to files in writeDirectory
#' @export
AnalyzeQuantiles <- function(readDirectory, writeDirectory, mean=TRUE, quantiles=NULL){

  #safety check
  if(!mean && is.null(quantiles)){
    stop("User needs to specify the mean or which quantiles to calculate. ")
  }

  #get files
  repFiles = list.dirs(path = readDirectory, full.names = TRUE, recursive = FALSE)
  patchFiles = lapply(X = repFiles, FUN = list.files, pattern = ".*\\.csv$")

  #subset females/males
  malePatches <- lapply(X = patchFiles, FUN = grep, pattern = "ADM", fixed = TRUE, value=TRUE)
  femalePatches <- lapply(X = patchFiles, FUN = grep, pattern = "AF1_Aggregate", fixed = TRUE, value=TRUE)

  #generate a list of all patches to run over
  patchList = unique(regmatches(x = patchFiles[[1]],
                                m = regexpr(pattern = "Patch[0-9]+",
                                            text = patchFiles[[1]],
                                            perl = TRUE)))

  #read in a file initially to get variables and setup return array
  testFile <- data.table::fread(input = file.path(repFiles[1], patchFiles[[1]][1]),
                                verbose = FALSE, showProgress = FALSE, drop = c("Time", "Patch"))


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



  #loop over all patches and do stats.
  for(patch in patchList){

    cat("Processing ", numReps, " reps of patch: ", patch,"\n",sep="")

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
                                                                 verbose = FALSE, showProgress = FALSE, drop = c("Time", "Patch")))
      popDataFemale[ ,repetition, ] <- as.matrix(data.table::fread(input = file.path(repFiles[repetition], femaleFiles[repetition]),
                                                                   verbose = FALSE, showProgress = FALSE, drop = c("Time", "Patch")))
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
                                file.path("ADM_Mean_", patch, ".csv", fsep = "")
      )
      femaleFileName <- file.path(writeDirectory,
                                  file.path("AF1_Mean_", patch, ".csv", fsep = "")
      )

      data.table::fwrite(x = as.data.frame(outputDataMale[ , ,1]),
                         file = maleFileName, col.names = TRUE, verbose = FALSE,
                         showProgress = FALSE, nThread = 1)
      data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,1]),
                         file = femaleFileName, col.names = TRUE, verbose = FALSE,
                         showProgress = FALSE, nThread = 1)
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
                                  file.path("ADM_Quantile_",
                                            formatC(x = quantiles[whichQuant], digits = 4,
                                                    format = "f", decimal.mark = "",
                                                    big.mark = NULL),
                                            "_", patch, ".csv", fsep = "")
        )
        femaleFileName <- file.path(writeDirectory,
                                    file.path("AF1_Quantile_",
                                              formatC(x = quantiles[whichQuant], digits = 4,
                                                      format = "f", decimal.mark = "",
                                                      big.mark = NULL),
                                              "_", patch, ".csv", fsep = "")
        )

        #write output
        data.table::fwrite(x = as.data.frame(outputDataMale[ , ,whichQuant]),
                           file = maleFileName, col.names = TRUE, verbose = FALSE,
                           showProgress = FALSE, nThread = 1)
        data.table::fwrite(x = as.data.frame(outputDataFemale[ , ,whichQuant]),
                           file = femaleFileName, col.names = TRUE, verbose = FALSE,
                           showProgress = FALSE, nThread = 1)
      }#end loop to write files
    }#end quantiles

    cat("Done with patch ", patch,".\n",sep="")

  }#end loop over patches
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
primePopMatrixArray <- function(primingMatrix,memoryWindow){
  array=vector(mode="list",length=memoryWindow)
  for(i in 1:length(array)){array[[i]]=primingMatrix}
  return(array)
}
