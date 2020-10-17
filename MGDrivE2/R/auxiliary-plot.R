################################################################################
#
#   MGDrivE2: Auxiliary plotting/data wrangling functions
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# Mosquitoes - Aquatic Stages
################################################################################
#######################################
# Genotype Designation
#######################################

#' Base Aquatic Function for Genotype Summary
#'
#' This function takes a given aquatic (egg, larval, pupal) stage and sums over
#' the Erlang-distributed stages, returning summary trajectories by genotype.
#'
#' This function is the base function for \code{\link{summarize_eggs_geno}},
#' \code{\link{summarize_larvae_geno}}, and \code{\link{summarize_pupae_geno}}.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{genotype}, and \code{value}.
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#' @param elp stage to summarize, one of: "egg", "larvae", "pupae"
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
base_aquatic_geno <- function(out,spn_P,elp){

  # get specific ELP status, remove nulls from human-only nodoes
  mosy_idx <- lapply(X = spn_P$ix, '[[', elp)
  mosy_idx[vapply(X = mosy_idx, FUN = is.null, FUN.VALUE = logical(length = 1))] <- NULL


  # get aux things
  # Since human-only nodes are possible, node labels aren't directly 1-to-n
  #  This uses the indices, grabs all of the column names, splits them into
  #  stage_geno_node, pulls out the nodes, gets the unique ones, then converts to integers
  # "if" protects 1-node case
  if(length(mosy_idx) == 1){
    nodes <- 1L
  } else {
    nodes <- as.integer(unique(vapply(X = strsplit(x = colnames(out)[unlist(mosy_idx)+1],
                                                   split = "_",fixed = TRUE),
                                      FUN = '[[', 3, FUN.VALUE = character(length = 1))))
  }

  num_nodes <- length(nodes)
  genos <- colnames(mosy_idx[[1]])
  numGeno <- length(genos)
  numRep <- dim(out)[3]



  # setup data holder
  mosyArray <- array(data = 0,dim = c(nrow(out),numGeno,num_nodes,numRep))


  # loop over reps
  for(r in 1:numRep){
    # loop over over nodes
    for(node in 1:num_nodes){
      # loop over genotypes
      for(gen in 1:numGeno){

        # consolodate over Erlang distributed stage duration
        mosyArray[ ,gen,node,r] <- rowSums(out[ ,mosy_idx[[node]][ ,gen] + 1,r])

      } # end geno loop
    } # end node loop
  } # end rep loop


  # setup return object
  retDF <- expand.grid("time" = out[ ,"time",1], "genotype" = genos,
                       "node" = nodes, "rep" = 1:numRep)

  # fill count data
  retDF$value <- as.vector(mosyArray)

  # if 1 node, remove node designation
  if(num_nodes == 1){retDF$node <- NULL}

  # if 1 rep, remove repetition designation
  if(numRep == 1){retDF$rep <- NULL}

  # return data frame
  return(retDF)
}


###################
# Genotype - Eggs
###################

#' Summarize Eggs by Genotype
#'
#' This function summarizes egg stage by genotype. It calls
#' \code{\link{base_aquatic_geno}} to do all of the work.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{genotype}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_eggs_geno <- function(out,spn_P){
  base_aquatic_geno(out = out, spn_P = spn_P, elp = "egg")
}


###################
# Genotype - Larvae
###################

#' Summarize Larvae by Genotype
#'
#' This function summarizes larval stage by genotype. It calls
#' \code{\link{base_aquatic_geno}} to do all of the work.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{genotype}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_larvae_geno <- function(out,spn_P){
  base_aquatic_geno(out = out, spn_P = spn_P, elp = "larvae")
}


###################
# Genotype - Pupae
###################

#' Summarize Pupal by Genotype
#'
#' This function summarizes pupal stage by genotype. It calls
#' \code{\link{base_aquatic_geno}} to do all of the work.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{genotype}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_pupae_geno <- function(out,spn_P){
  base_aquatic_geno(out = out, spn_P = spn_P, elp = "pupae")
}


#######################################
# Erlang-Stage Designation
#######################################

#' Base Aquatic Function for Erlang-Stage Summary
#'
#' This function takes a given aquatic (egg, larval, pupal) stage and sums over
#' the genotypes, returning summary trajectories by Erlang-distributed stage.
#'
#' This function is the base function for \code{\link{summarize_eggs_stage}},
#' \code{\link{summarize_larvae_stage}}, and \code{\link{summarize_pupae_stage}}.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{Erlang-stage}, and \code{value}.
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#' @param elp stage to summarize, one of: "egg", "larvae", "pupae"
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
base_aquatic_stage <- function(out,spn_P,elp){

  # get specific ELP status, remove nulls from human-only nodoes
  mosy_idx <- lapply(X = spn_P$ix, '[[', elp)
  mosy_idx[vapply(X = mosy_idx, FUN = is.null, FUN.VALUE = logical(length = 1))] <- NULL


  # get aux things
  # Since human-only nodes are possible, node labels aren't directly 1-to-n
  #  This uses the indices, grabs all of the column names, splits them into
  #  stage_geno_node, pulls out the nodes, gets the unique ones, then converts to integers
  # "if" protects 1-node case
  if(length(mosy_idx) == 1){
    nodes <- 1L
  } else {
    nodes <- as.integer(unique(vapply(X = strsplit(x = colnames(out)[unlist(mosy_idx)+1],
                                                   split = "_",fixed = TRUE),
                                      FUN = '[[', 3, FUN.VALUE = character(length = 1))))
  }

  num_nodes <- length(nodes)
  erlangStage <- rownames(mosy_idx[[1]])
  numES <- length(erlangStage)
  numRep <- dim(out)[3]



  # setup data holder
  mosyArray <- array(data = 0,dim = c(nrow(out),numES,num_nodes,numRep))


  # loop over reps
  for(r in 1:numRep){
    # loop over over nodes
    for(node in 1:num_nodes){
      # loop over genotypes
      for(stage in 1:numES){

        # consolodate over Erlang distributed stage duration
        mosyArray[ ,stage,node,r] <- rowSums(out[ ,mosy_idx[[node]][stage, ] + 1,r])

      } # end geno loop
    } # end node loop
  } # end rep loop


  # setup return object
  retDF <- expand.grid("time" = out[ ,"time",1], "Erlang-stage" = erlangStage,
                       "node" = nodes, "rep" = 1:numRep)

  # fill count data
  retDF$value <- as.vector(mosyArray)

  # if 1 node, remove node designation
  if(num_nodes == 1){retDF$node <- NULL}

  # if 1 rep, remove repetition designation
  if(numRep == 1){retDF$rep <- NULL}

  # return data frame
  return(retDF)
}


###################
# Erlang-Stage - Eggs
###################

#' Summarize Eggs by Erlang-Stage
#'
#' This function summarizes egg stage by Erlang-stages. It calls
#' \code{\link{base_aquatic_stage}} to do all of the work.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{Erlang-stage}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_eggs_stage <- function(out,spn_P){
  base_aquatic_stage(out = out, spn_P = spn_P, elp = "egg")
}


###################
# Erlang-Stage - Larvae
###################

#' Summarize Larval by Erlang-Stage
#'
#' This function summarizes larval stage by Erlang-stages. It calls
#' \code{\link{base_aquatic_stage}} to do all of the work.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{Erlang-stage}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_larvae_stage <- function(out,spn_P){
  base_aquatic_stage(out = out, spn_P = spn_P, elp = "larvae")
}


###################
# Erlang-Stage - Pupae
###################

#' Summarize Pupal by Erlang-Stage
#'
#' This function summarizes pupal stage by Erlang-stages. It calls
#' \code{\link{base_aquatic_stage}} to do all of the work.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{Erlang-stage}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_pupae_stage <- function(out,spn_P){
  base_aquatic_stage(out = out, spn_P = spn_P, elp = "pupae")
}


################################################################################
# Mosquitoes - Female
################################################################################
#######################################
# Lifecycle
#######################################

#' Summarize Adult Females (One Node or Metapopulation Network, Lifecycle Model)
#'
#' For \code{MGDrivE2} simulations of mosquito lifecycle dynamics in a single node
#' or metapopulation network, this function sums over the male mate genotype to
#' get population trajectories of adult female mosquitoes by their genotype.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}} or \code{\link{spn_P_lifecycle_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{genotype}, and \code{value}.
#'
#' For examples of using this function, this or any vignette which visualizes output:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_females <- function(out,spn_P){

  # get constants for later
  gNames <- dimnames(spn_P$ix[[1]]$females)[[1]]
  numGeno <- dim(spn_P$ix[[1]]$females)[1]
  num_nodes <- length(spn_P$ix)
  numRep <- dim(out)[3]



  # setup data holder
  fArray <- array(data = 0, dim = c(nrow(out),numGeno,num_nodes,numRep))

  # loop over reps, nodes/genotypes, collapse females by mate
  for(r in 1:numRep){
    for(node in 1:num_nodes){
      for(gen in 1:numGeno){
        fArray[ ,gen,node,r] <- rowSums(out[,spn_P$ix[[node]]$females[gen, ,1]+1,r])
      }
    }
  }


  # setup return df
  retDF <- expand.grid("time" = out[ ,"time",1], "genotype" = gNames,
                       "node" = 1:num_nodes, "rep" = 1:numRep)

  # fill count data
  retDF$value <- as.vector(fArray)

  # if 1 node, remove node designation
  if(num_nodes == 1){retDF$node <- NULL}

  # if 1 rep, remove repetition designation
  if(numRep == 1){retDF$rep <- NULL}

  # return data frame
  return(retDF)
}


#######################################
# Epi - SEI
#######################################

#' Summarize Adult Females (One Node or Metapopulation Network, SEI Mosquitoes)
#'
#' For \code{MGDrivE2} simulations of mosquito epidemiological dynamics in a single
#' node or metapopulation network, this function sums over the male mate genotype
#' as well as EIP bins to get population trajectories of adult female mosquitoes
#' by their genotype and (S,E,I) status.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{inf}, \code{genotype}, and \code{value}.
#'
#' For examples of using this function, this or any vignette which simulates epi dynamics:
#' \code{vignette("epi-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param spn_P the places of the SPN, see details
#'
#' @return a 4 to 6 column dataframe for plotting with ggplot2
#'
#' @export
summarize_females_epi <- function(out,spn_P){

  # get females for epi network, then remove nulls
  f_idx <- lapply(X = spn_P$ix, '[[', 'females')
  f_idx[vapply(X = f_idx, FUN = is.null, FUN.VALUE = logical(length = 1))] <- NULL

  # get aux things
  # Since human-only nodes are possible, node labels aren't directly 1-to-n
  #  This uses the indices, grabs all of the column names, splits them into
  #  F_geno_geno_stage_node, pulls out the nodes, gets the unique ones, then converts to integers
  # "if" protects 1-node case
  if(length(f_idx) == 1){
    nodes <- 1L
  } else {
    nodes <- as.integer(unique(vapply(X = strsplit(x = colnames(out)[unlist(f_idx)+1],
                                                   split = "_",fixed = TRUE),
                                      FUN = '[[', 5, FUN.VALUE = character(length = 1))))
  }

  num_nodes <- length(nodes)
  genos <- dimnames(f_idx[[1]])[[1]]
  numGeno <- length(genos)
  numRep <- dim(out)[3]

  statusSEI <- dimnames(f_idx[[1]])[[3]]

  expCols <- 2:(length(statusSEI)-1)

  # setup data holder
  #  this only works for SEI mosquitoes, because that gives "3" states here
  fArray <- array(data = 0,dim = c(nrow(out),numGeno * 3,num_nodes,numRep))


  # loop over reps
  for(r in 1:numRep){
    # loop over over nodes
    for(node in 1:num_nodes){
      # loop over genotypes
      genStep <- 0
      for(gen in 1:numGeno){
        # consolodate 3 things
        # susceptible
        fArray[ ,genStep + 1,node,r] <- rowSums(out[ ,f_idx[[node]][gen, ,"S"] + 1,r])
        # exposed
        fArray[ ,genStep + 2,node,r] <- rowSums(out[ ,f_idx[[node]][gen, ,expCols] + 1,r])
        # infectious
        fArray[ ,genStep + 3,node,r] <- rowSums(out[ ,f_idx[[node]][gen, ,"I"] + 1,r])
        # increment step
        genStep <- genStep + 3

      } # end geno loop
    } # end node loop
  } # end rep loop


  # setup return object
  retDF <- expand.grid("time" = out[ ,"time",1], "inf" = c("S","E","I"),
                       "genotype" = genos, "node" = nodes, "rep" = 1:numRep)

  # fill count data
  retDF$value <- as.vector(fArray)

  # if 1 node, remove node designation
  if(num_nodes == 1){retDF$node <- NULL}

  # if 1 rep, remove repetition designation
  if(numRep == 1){retDF$rep <- NULL}

  # return data frame
  return(retDF)
}


################################################################################
# Mosquitoes - Male
################################################################################

#' Summarize Adult Males (One Node or Metapopulation Network)
#'
#' For \code{MGDrivE2} simulations of mosquito lifecycle dynamics or human infection
#' dynamics, in a node or metapopulation network, this function summarizes
#' population trajectories of adult male mosquitoes by their genotype.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{genotype}, and \code{value}.
#'
#' For examples of using this function, this or any vignette which visualizes output:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#'
#' @return a 3 to 5 column dataframe for plotting with ggplot2
#'
#' @export
summarize_males <- function(out){

  # get names of all states
  states <- colnames(out)
  numRep <- dim(out)[3]


  # grab only the male mosquito states
  #  identifier is M, anchored at the front, avoiding genos and node
  m_idx <- grep(pattern = "^M_", x = states)

  # split states into M, genotypes, and node
  gens_nodes <- strsplit(x = states[m_idx], split = "_", fixed = TRUE)

  # safety net if one node
  if(all(lengths(gens_nodes) < 3)){
    nodes = 1
  } else {
    nodes <- as.integer(x = vapply(X = gens_nodes, '[[',3, FUN.VALUE = character(length = 1)))
  }

  # get genotypes
  genos <- vapply(X = gens_nodes, '[[',2, FUN.VALUE = character(length = 1))

  # setup return df
  retDF <- expand.grid("time" = out[ ,"time",1], "genotype" = unique(genos),
                       "node" = unique(nodes), "rep" = 1:numRep)

  # fill count data
  retDF$value <- as.vector(out[ ,m_idx, ])

  # if 1 node, remove node designation
  if(length(nodes) == 1){retDF$node <- NULL}

  # if 1 rep, remove repetition designation
  if(numRep == 1){retDF$rep <- NULL}

  # return data frame
  return(retDF)
}


################################################################################
# Humans
################################################################################

#' Base Function for Human Summary
#'
#' This function takes a given infection ('S','E','I','R') status and returns a
#' summary trajectory
#'
#' This function is the base function for \code{\link{summarize_humans_epiSIS}},
#' \code{\link{summarize_humans_epiSEIR}}.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{inf}, \code{genotype}, and \code{value}.
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#' @param infState type of humans to summarize: 'S','E','I','R'
#'
#' @return a 4 to 6 column dataframe for plotting with ggplot2
#'
base_summarize_humans <- function(out,infState){

  # get names of all states
  states <- colnames(out)
  numRep <- dim(out)[3]

  # grab only the human states
  h_idx <- grep(pattern = paste0(c("^H_[",infState,"](_[0-9]+)?"),collapse = ""),
               x = states)

  # split states into infection state, human, and node
  inf_nodes <- strsplit(x = states[h_idx], split = "_", fixed = TRUE)

  # safety net if one node
  if(all(lengths(inf_nodes) < 3)){
    nodes = 1
  } else {
    nodes <- as.integer(x = vapply(X = inf_nodes, '[[',3, FUN.VALUE = character(length = 1)))
  }

  # setup return df
  retDF <- expand.grid("time" = out[ ,"time",1], "inf" = infState, "genotype" = "human",
                        "node" = unique(nodes), "rep" = 1:numRep)

  # fill count data
  retDF$value <- as.vector(out[ ,h_idx, ])

  # if 1 node, remove node designation
  if(length(nodes) == 1){retDF$node <- NULL}

  # if 1 rep, remove repetition designation
  if(numRep == 1){retDF$rep <- NULL}

  # return data frame
  return(retDF)
}


#######################################
# Humans - SIS
#######################################

#' Summarize Humans (One Node or Metapopulation Network, SEI Mosquitoes - SIS Humans)
#'
#' For \code{MGDrivE2} simulations of mosquito epidemiological dynamics in a
#' node or network, this function summarizes human infection status, S and I. It
#' uses \code{\link{base_summarize_humans}} to do all of the work.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{inf}, \code{genotype}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("epi-node", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#'
#' @return a 4 to 6 column dataframe for plotting with ggplot2
#'
#' @export
summarize_humans_epiSIS <- function(out){
  base_summarize_humans(out = out, infState = c("S","I"))
}


#######################################
# Humans - SEIR
#######################################

#' Summarize Humans (One Node or Metapopulation Network, SEI Mosquitoes - SEIR Humans)
#'
#' For \code{MGDrivE2} simulations of mosquito epidemiological dynamics in a
#' node or network, this function summarizes human infection status, S, E, I, and R.
#' It uses \code{\link{base_summarize_humans}} to do all of the work.
#'
#' The return object depends on the data provided. If the simulation was only 1 node,
#' then no \code{node} designation is returned. If only one repetition was performed,
#' no \code{rep} designation is returned. Columns always returned include: \code{time},
#' \code{inf}, \code{genotype}, and \code{value}.
#'
#' For examples of using this function, see:
#' \code{vignette("seir-dynamics", package = "MGDrivE2")}
#'
#' @param out the output of \code{\link[MGDrivE2]{sim_trajectory_R}}
#'
#' @return a 4 to 6 column dataframe for plotting with ggplot2
#'
#' @export
summarize_humans_epiSEIR <- function(out){
  base_summarize_humans(out = out, infState = c("S","E","I","R"))
}
