###############################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   MGDrivE Graphics
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################
#######################################
# Color Utility
#######################################

#' Utility to Imitate ggplot2 Colors
#'
#' Sample at equally spaced intervals along the color wheel
#'
#' @param n Number of colors
#' @param alpha Transparency
#'
#' @importFrom grDevices hcl
#'
ggColUtility <- function(n, alpha = .75) {
  # equally spaced intervals
  hues = seq(15, 375, length = n + 1)

  # return colors
  return(hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n])
}

#######################################
# Plot 1 Trace
#######################################

#' Plot
#'
#' Plots one run from MGDrivE
#'
#' @usage plotMGDrivESingle(readDir, whichPatches = NULL, totalPop = FALSE,
#'                          nonZeroGen = FALSE, lwd = 0.75, alpha = 0.75)
#'
#' @param readDir Path to file from single-run of MGDrivE or from analysis function
#' @param whichPatches Vector of patches to plot, must be less than 15. Default is NULL if less than 15 patches
#' @param totalPop Boolean, to plot the total population or not.
#' @param nonZeroGen Boolean, to plot genotypes that are always zero in simulation
#' @param lwd Double, specify the line width for plotting
#' @param alpha Double, specify the opacity for plotting
#'
#' @importFrom graphics box grid layout legend matplot mtext par plot.new title
#' @importFrom utils tail
#'
#' @details This function plots output from one run or one set of runs after being
#' analyzed. Setting totalPop to FALSE keeps it from plotting the total population.
#' NonZeroGen accounts for genotypes that could exist, but are not created in the
#' simulation. Default is FALSE, as this is easier to read on a plot.
#'
#' @examples
#' \dontrun{
#' # Requires the user to have run MGDrivE, deterministic or stochastic, analyzed
#' #  the data, and stored it in the directory shown below.
#' # See vignette for complete example
#'
#' # Folder where single run is stored
#' fPath <- "path/to/data/containing/folder"
#'
#' # plot output to see effect
#' plotMGDrivESingle(readDir=fPath,totalPop = TRUE,lwd=3.5,alpha=1)
#' }
#'
#' @export
plotMGDrivESingle <- function(readDir, whichPatches = NULL, totalPop = FALSE,
                              nonZeroGen = FALSE, lwd = 0.75, alpha = 0.75){

  #keep old plot parameters to reset later
  oldPar <- par(no.readonly = TRUE)
  on.exit(expr = par(oldPar)) #reset par()

  ####################
  # get files to plot, check number of patches
  ####################
  mFiles = list.files(path = readDir, pattern = "^M_.*csv$", full.names = TRUE)
  fFiles = list.files(path = readDir, pattern = "^F_.*csv$", full.names = TRUE)

  # Test that number of patches matches what we can plot
  patches = unique(x = regmatches(x = mFiles, m = regexpr(pattern = "Patch[0-9]+", text = mFiles)))

  # check if user chose specific patches
  if(length(patches)>15 && is.null(whichPatches)){
    stop("There are more than 15 patches in the simulation.
         Please select less than 15 to plot.")
  }
  if(!is.null(whichPatches)){
    # make sure not too many. I dont' know what that number is, but at some point
    #  the plotting function breaks down.
    if(length(whichPatches)>15){
      stop("Please select less than 15 patches.")
    }
    # select the patches, if they select less than the total number of patches.
    if(length(whichPatches)<=length(patches)){
      patches = patches[whichPatches]
    }
  }

  numPatches <- length(patches)


  ####################
  # scan test file for names and existing genotypes
  ####################
  # test file to get sizes
  # THIS ASSUMES MALE/FEMALE FILE ORDER. MAY CHANGE IF WE RENAME FILES
  columnNames <- scan(file = mFiles[1], what = character(), sep = ",", quiet = TRUE, nlines = 1)

  mFile <- matrix(data = scan(file = mFiles[1], what = double(), sep = ",", skip = 1, quiet = TRUE),
                     ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))

  fFile <- matrix(data = scan(file = fFiles[1], what = double(), sep = ",", skip = 1, quiet = TRUE),
                  ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))


  # base genotypes
  genotypes <- columnNames
  # if non-zero gens only
  #   test male and female files, in case of sex-specific drive
  #   Does only test first patch, which if releases aren't done there, could be wrong.
  if(!nonZeroGen){genotypes <- genotypes[colSums(mFile[ ,genotypes])!=0 | colSums(fFile[ ,genotypes])!=0]}


  ####################
  # Read in all data, subset by desired genotypes
  ####################
  # create data list holder, and then fill it
  maleData=femaleData=array(data = 0L, dim = c(nrow(mFile),length(genotypes)+1,numPatches),
                            dimnames = list(NULL, c(genotypes, "Total"), patches))

  for(patch in patches){
    # get male file name
    names = grep(pattern = patch, x = mFiles, fixed = TRUE, value = TRUE)

    mFile[] = matrix(data = scan(file = names, what = double(),
                                 sep = ",", skip = 1, quiet = TRUE),
                     ncol = length(columnNames), byrow = TRUE)
    maleData[,,patch] = cbind(mFile[ ,genotypes], rowSums(x = mFile[ ,-1]))

    # get female file name
    names = grep(pattern = patch, x = fFiles, fixed = TRUE, value = TRUE)

    fFile[] = matrix(data = scan(file = names, what = double(),
                                 sep = ",", skip = 1, quiet = TRUE),
                     ncol = length(columnNames), byrow = TRUE)
    femaleData[,,patch] = cbind(fFile[ ,genotypes], rowSums(x = fFile[ ,-1]))

  }# end reading in data

  ####################
  # setup colors and final genotype size
  ####################
  # test for total pop
  if(totalPop){
    genotypes <- c(genotypes, "Total")
  }

  # remove time
  genotypes <- genotypes[-1]

  # get color swath
  col <- ggColUtility(n = length(genotypes), alpha=alpha)


  ####################
  # plot layout
  ####################
  lmatrix <- matrix(data = 1:(numPatches*3), nrow = numPatches, ncol = 3, byrow = TRUE)
  if(numPatches>1){
    # fill in rest of plot labels
    lmatrix[2:numPatches, c(1,2)] <- matrix(data = 4:(3+2*(numPatches-1)),
                                            ncol = 2, byrow = TRUE)
    # legend gets whole right side
    lmatrix[,3] <- 3
  }

  layout(lmatrix, widths = c(3,3,1))

  ####################
  # plot!
  ####################
  # plot first patch and the legend
  # female
  par(mar = c(2,3,3,1), las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 1.75)
  matplot(x = femaleData[ ,1,patches[1],drop = FALSE], y = femaleData[ , genotypes, patches[1]],
          type = "l", lty = 1, main = "Females", ylab = "", lwd=lwd,
          ylim = c(0, 1.1*max(femaleData[ , genotypes, patches[1]])),
          xlim = c(0, tail(x = femaleData[ ,1,1], n = 1)), yaxs = "i", xaxs = "i",
          col = col, panel.first=grid())
  title(ylab = "Population\n", line = 2)
  box(lwd = 2)

  # male
  par(mar = c(2,2,3,1), las = 1)
  matplot(x = maleData[ ,1, patches[1],drop = FALSE], y = maleData[ , genotypes, patches[1]],
          type = "l", lty = 1, main = "Males", ylab = "", lwd=lwd,
          ylim = c(0, 1.1*max(maleData[ , genotypes, patches[1]])),
          xlim = c(0, tail(x = maleData[ ,1,1], n = 1)), yaxs = "i", xaxs = "i",
          col = col, panel.first=grid())
  mtext(patches[1], side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
  box(lwd = 2)

  # legend
  plot.new()
  par(font = 2)
  legend(x = 'left',
         inset = 0, # will have to play with this as widths change
         seg.len = 2, # length of line denoting the colors
         x.intersp = 0.9,
         y.intersp = 0.9, # vertical space between lines
         title = 'Genotypes',
         legend = genotypes,
         col = col,
         bty = "n",
         #bg = "lightblue",
         xpd = TRUE,
         pch = 15,
         pt.cex = 2,
         # lty = c(1,1,1,1),
         # lwd=16,
         cex = 1)

  # rest of the patches
  if(numPatches>1){
    for(patch in patches[-1]){
      # female
      par(mar = c(2,3,1,1), las = 1)
      matplot(x = femaleData[ ,1, patch, drop = FALSE], y = femaleData[ , genotypes, patch],
              type = "l", lty = 1, ylab = "", lwd=lwd,
              ylim = c(0, 1.1*max(femaleData[ , genotypes, patch])),
              xlim = c(0, tail(x = femaleData[ ,1,1], n = 1)), yaxs = "i", xaxs = "i",
              col = col, panel.first=grid())
      title(ylab = "Population\n", line = 2)
      box(lwd = 2)

      # male
      par(mar = c(2,2,1,1))
      matplot(x = maleData[ ,1, patch, drop = FALSE], y = maleData[ , genotypes, patch],
              type = "l", lty = 1, ylab = "", lwd=lwd,
              ylim = c(0, 1.1*max(maleData[ , genotypes, patch])),
              xlim = c(0, tail(x = maleData[ ,1,1], n = 1)), yaxs = "i", xaxs = "i",
              col = col, panel.first=grid())
      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
    } # end patch loop
  } # end rest of patches

}

#######################################
# Plot multiple Trace
#######################################

#' Plot
#'
#' Plots several traces from MGDrivE, assuming each set is another repetition
#' from the same experiment. \cr
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
#'    * ...
#'
#' @usage plotMGDrivEMult(readDir, whichPatches = NULL, totalPop = FALSE,
#'                        nonZeroGen = FALSE, lwd = 0.75, alpha = 0.75)
#'
#' @param readDir Directory to find repetition folders in
#' @param whichPatches Vector of patches to plot, must be less than 15. Default is NULL if less than 15 patches
#' @param totalPop Boolean, to plot the total population or not. Default is FALSE
#' @param nonZeroGen Boolean, to plot genotypes that are always zero in simulation
#' @param lwd Double, specify the line width for plotting
#' @param alpha Double, specify the opacity for plotting
#'
#' @importFrom graphics matlines
#'
#' @details This function plots output from one run or one set of runs after
#' being analyzed. Setting totalPop to FALSE keeps it from plotting the total
#' population. NonZeroGen accounts for genotypes that could exist, but are not
#' created in the simulation. Default is FALSE, as this is easier to read on a plot.
#'
#' @examples
#' \dontrun{
#' # Requires the user to have run MGDrivE, logically stochastic, analyzed
#' #  the data, and stored it in the directory shown below.
#' # See vignette for complete example
#'
#' # Folder where single run is stored
#' fPath <- "path/to/data/containing/folder"
#'
#' # plot output to see effect
#' plotMGDrivEMult(readDir=fPath,totalPop = TRUE,lwd=3.5,alpha=1)
#' }
#'
#' @export
plotMGDrivEMult <- function(readDir, whichPatches = NULL, totalPop = FALSE,
                            nonZeroGen = FALSE, lwd = 0.75, alpha = 0.75){

  #keep old plot parameters to reset later
  oldPar <- par(no.readonly = TRUE)
  on.exit(expr = par(oldPar)) #reset par()

  ####################
  # check number of patches
  ####################
  mFiles = lapply(X = list.dirs(path = readDir, recursive = FALSE),
                    FUN = list.files, pattern = "^M_.*csv$", full.names = TRUE)
  fFiles = lapply(X = list.dirs(path = readDir, recursive = FALSE),
                  FUN = list.files, pattern = "^F_.*csv$", full.names = TRUE)

  # Test that number of patches matches what we can plot
  patches = unique(x = regmatches(x = mFiles[[1]], m = regexpr(pattern = "Patch[0-9]+", text = mFiles[[1]])))

  # check if user chose specific patches
  if(length(patches)>15 && is.null(whichPatches)){
    stop("There are more than 15 patches in the simulation.
         Please select less than 15 to plot.")
  }
  if(!is.null(whichPatches)){
    # make sure not too many. I dont' know what that number is, but at some point
    #  the plotting function breaks down.
    if(length(whichPatches)>15){
      stop("Please select less than 15 patches.")
    }
    # select the patches, if they select less than the total number of patches.
    if(length(whichPatches)<=length(patches)){
      patches = patches[whichPatches]
    }
  }

  numPatches <- length(patches)
  numReps <- length(mFiles)

  ####################
  # scan test file for names and existing genotypes
  ####################
  # test file to get sizes
  columnNames <- scan(file = mFiles[[1]][1], what = character(), sep = ",", quiet = TRUE, nlines = 1)

  mFile <- matrix(data = scan(file = mFiles[[1]][1], what = double(), sep = ",", skip = 1, quiet = TRUE),
                  ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))

  fFile <- matrix(data = scan(file = fFiles[[1]][1], what = double(), sep = ",", skip = 1, quiet = TRUE),
                  ncol = length(columnNames), byrow = TRUE, dimnames = list(NULL,columnNames))


  # base genotypes
  genotypes <- columnNames
  # if non-zero gens only
  #   test male and female files, in case of sex-specific drive
  #   Does only test first patch, which if releases aren't done there, could be wrong.
  if(!nonZeroGen){genotypes <- genotypes[colSums(mFile[ ,genotypes])!=0 | colSums(fFile[ ,genotypes])!=0]}

  ####################
  # Read in all data, subset by desired genotypes
  ####################
  # create data list holder, and then fill it
  #  These holders store every patch for every repetition.
  #  The list length is the number of repetitions, the array inside is each patch
  maleData=femaleData=rep(x = list(array(data = 0L, dim = c(nrow(mFile),length(genotypes)+1,numPatches),
                                         dimnames = list(NULL, c(genotypes, "Total"), patches))),
                          numReps)

  # loop over repetition
  for(nRep in 1:numReps){
    # loop over patch
    for(patch in patches){
      # get male patch name
      names = grep(pattern = patch, x = mFiles[[nRep]], fixed = TRUE, value = TRUE)

      # assumes male is before female alphabetically
      mFile[] = matrix(data = scan(file = names, what = double(),
                                   sep = ",", skip = 1, quiet = TRUE),
                       ncol = length(columnNames), byrow = TRUE)
      maleData[[nRep]][,,patch] = cbind(mFile[ ,genotypes], rowSums(x = mFile[ ,-1]))


      # get female patch name
      names = grep(pattern = patch, x = fFiles[[nRep]], fixed = TRUE, value = TRUE)

      fFile[] = matrix(data = scan(file = names, what = double(),
                                   sep = ",", skip = 1, quiet = TRUE),
                       ncol = length(columnNames), byrow = TRUE)
      femaleData[[nRep]][,,patch] = cbind(fFile[ ,genotypes], rowSums(x = fFile[ ,-1]))

    }# end patches
  }# end repetitions

  ####################
  # setup colors and final genotype size
  ####################
  # test for total pop
  if(totalPop){
    genotypes <- c(genotypes, "Total")
  }

  # remove time
  genotypes <- genotypes[-1]

  # get color swath
  col <- ggColUtility(n = length(genotypes), alpha=alpha)

  ####################
  # plot layout
  ####################
  lmatrix <- matrix(data = 1:(numPatches*3), nrow = numPatches, ncol = 3, byrow = TRUE)
  if(numPatches>1){
    #fill in rest of plot labels
    lmatrix[2:numPatches, c(1,2)] <- matrix(data = 4:(3+2*(numPatches-1)),
                                            ncol = 2, byrow = TRUE)
    #legend gets whole right side
    lmatrix[,3] <- 3
  }

  layout(lmatrix, widths = c(3,3,1))

  ####################
  # plot!
  ####################
  # plot first patch and the legend
  #  female
  #  set parameters
  par(mar = c(2,3,3,1), las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 1.75)
  # plot first repetition
  matplot(x = femaleData[[1]][ , 1, patches[1], drop = FALSE], y = femaleData[[1]][ , genotypes, patches[1]],
          type = "l", lty = 1, main = "Females", ylab = "", lwd=lwd,
          ylim = c(0, 1.1*max(femaleData[[1]][ , genotypes, patches[1]])),
          xlim = c(0, tail(x = femaleData[[1]][ , 1, 1], n = 1)), yaxs = "i", xaxs = "i",
          col = col, panel.first=grid())
  # add remaining repetitions to the plot
  for(nRep in 2:numReps){
    matlines(x = femaleData[[nRep]][ ,1, patches[1]], y = femaleData[[nRep]][ , genotypes, patches[1]],
             type = "l", lty = 1,lwd=lwd, col = col)
  }
  title(ylab = "Population\n", line = 2)
  box(lwd = 2)

  # Male
  # set parameters
  par(mar = c(2,2,3,1), las = 1)
  # plot first repetition
  matplot(x = maleData[[1]][ ,1, patches[1], drop = FALSE], y = maleData[[1]][ , genotypes, patches[1]],
          type = "l", lty = 1, main = "Males", ylab = "", lwd=lwd,
          ylim = c(0, 1.1*max(maleData[[1]][ , genotypes, patches[1]])),
          xlim = c(0, tail(x = maleData[[1]][ , 1, 1], n = 1)), yaxs = "i", xaxs = "i",
          col = col, panel.first=grid())
  # add remaining repetitions to the plot
  for(nRep in 2:numReps){
    matlines(x = maleData[[nRep]][ ,1, patches[1], drop = FALSE], y = maleData[[nRep]][ , genotypes, patches[1]],
             type = "l", lty = 1,lwd=lwd, col = col)
  }
  mtext(patches[1], side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
  box(lwd = 2)

  # legend
  plot.new()
  par(font = 2)
  legend(x = 'left',
         inset = 0, # will have to play with this as widths change
         seg.len = 2, # length of line denoting the colors
         x.intersp = 0.9,
         y.intersp = 0.9, # vertical space between lines
         title = 'Genotypes',
         legend = genotypes,
         col = col,
         bty = "n",
         # bg = "lightblue",
         xpd = TRUE,
         pch = 15,
         pt.cex = 2,
         # lty = c(1,1,1,1),
         # lwd=16,
         cex = 1)

  # rest of the patches
  if(numPatches>1){
    # loop over patches
    for(patch in patches[-1]){
      # female
      par(mar = c(2,3,1,1), las = 1)
      matplot(x = femaleData[[1]][ ,1, patch], y = femaleData[[1]][ , genotypes, patch],
              type = "l", lty = 1, ylab = "", lwd=lwd,
              ylim = c(0, 1.1*max(femaleData[[1]][ , genotypes, patch])),
              xlim = c(0, tail(x = femaleData[[1]][ , 1, patch], n = 1)), yaxs = "i", xaxs = "i",
              col = col, panel.first=grid())
      for(nRep in 2:numReps){
        matlines(x = femaleData[[nRep]][ ,1, patch], y = femaleData[[nRep]][ , genotypes, patch],
                 type = "l", lty = 1,lwd=lwd, col = col)
      }
      title(ylab = "Population\n", line = 2)
      box(lwd = 2)

      # male
      par(mar = c(2,2,1,1))
      matplot(x = maleData[[1]][ ,1, patch], y = maleData[[1]][ , genotypes, patch],
              type = "l", lty = 1, ylab = "", lwd=lwd,
              ylim = c(0, 1.1*max(maleData[[1]][ , genotypes, patch])),
              xlim = c(0, tail(x = maleData[[1]][ , 1, patch], n = 1)), yaxs = "i", xaxs = "i",
              col = col, panel.first=grid())
      for(nRep in 2:numReps){
        matlines(x = maleData[[nRep]][ ,1, patch], y = maleData[[nRep]][ , genotypes, patch],
                 type = "l", lty = 1,lwd=lwd, col = col)
      }

      mtext(patch, side = 4, line = 0.5, las = 0, cex = 0.9, font = 2)
      box(lwd = 2)
    } # end patch loop
  } # end rest of patches

}
