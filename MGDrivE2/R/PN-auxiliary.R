################################################################################
#
#   MGDrivE2: internal functions to make epi metapopulation network models
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2019
#
################################################################################


################################################################################
#   stitch together transitions {T} for the metapop:
#   helper functions
################################################################################

# checker functions for movement, vectorized
h_move_check <- function(x,y){ (x %in% c("b","h")) & (y %in% c("b","h")) }

m_move_check <- function(x,y){ (x %in% c("b","m")) & (y %in% c("b","m")) }

# function to check that user-input movement is legal
check_move_legal <- function(node_list,move,checker,hm){
  # make sure right checker is called
  if(!(hm %in% c("human","mosquito"))){stop("check_move_legal: bad input for arg 'hm'")}

  # check combos, same shape/order as move
  move_legal <- outer(X = node_list,Y = node_list,FUN = checker)

  # get indices of non-viable moves
  no_move <- which(x = !move_legal, arr.ind = TRUE, useNames = FALSE)

  if(any(xor(move[no_move],move_legal[no_move]))){
    stop(paste0("error: movement matrix entry is not allowed, please fix ",
                hm," movement matrix") )
  }

}
