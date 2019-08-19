#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************
 * Auxilary functions to get rid of other library dependencies
 ******************************************************************************/

// //' Make a Symmetric Cube
// //'
// //' This function makes a lower-triangular cube symmetric over the z-axis.
// //' It was written to remove the dependency "Matrix".
// //' It is only used in building cubes.
// //'
// //' @param lowerMat A lower-triangular matrix of depth 1 or more
// //'
// // [[Rcpp::export]]
// void symCube(arma::dcube& lowerMat){
//
//   for(int i=0; i<lowerMat.n_slices; i++){
//     lowerMat.slice(i) = symmatl(lowerMat.slice(i));
//   }
//
// }

//' Shift a Vector
//'
//' Shift a population vector by one day and insert the new population.
//' This was written to remove the dependency "binhf".
//'
//' @param popVector List of population vectors of length(Tegg+Tlarva+Tpupa)
//' @param newPop Vector of length equal to the number of genotypes
//'
// [[Rcpp::export]]
void shiftAndUpdatePopVector( ListOf<NumericVector>& popVector, const NumericVector newPop){
  //This  function shifts the list of population vectors one to the right, then
  // updates the first position with the new population.
  // Only shifts things by 1 point!
  // works on doubles and ints
  // works on lists length 1
  // OVERWRITES THE INPUT!!! BE CAREFUL

  /*This function is used in Patch-Methods oneDay_UpdatePopulation_AgeClass()
   *Originally, it was meant to replace a call to shift() from the binhf package.
   *It became easier to remake the whole function. It takes the population list,
   * shifts it to the right by 1, then updates the first slot with the new
   * population. It won't move things to the left, or by a step other than one,
   * and the safety checks are gone. FYI.
   */

  //Shift to the right, starting from the second to last element
  for(int i = popVector.size()-2; i>=0; i--){

    popVector[i+1] = popVector[i];

  }

  //set the first position with the new population
  popVector[0] = newPop;
}

//' Dirichlet Distribution
//'
//' Make a single draw from a Dirichlet distribution with the shape parameter
//' one. This replaces the MCMCpack rDirichlet function, which was wholly written
//' in R.
//'
//' @param migrationPoint Vector of weights for draws. Must be positive.
//'
// [[Rcpp::export]]
NumericVector rDirichlet(const NumericVector& migrationPoint){

  //set up return things
  NumericVector probs(migrationPoint.length());

  //This is a Dirichlet distribtuion
  for(int i = 0; i<migrationPoint.length(); i++){
    probs[i] = R::rgamma(migrationPoint[i], 1.0);
  }
  probs = probs / sum(probs);

  //return
  return(probs);
}

//' Quantiles Function
//'
//' Calculate the given quantiles of a matrix.
//'
//' @usage quantileC(Trials, Probs)
//'
//' @param Trials Integer matrix to calculate quantiles over
//' @param Probs Vector of quantiles
//'
//' @details This function calculates the given quantiles over the rows of an
//' integer matrix. It uses method 8 of the stat::quantiles() function. It gives
//' the same result, to numerical accuracy, and is designed to handle matrix input.
//' It is only designed to work on integer matrices!
//'
//' @return Numeric Matrix
//'
// [[Rcpp::export]]
NumericMatrix quantileC(IntegerMatrix& Trials, const NumericVector& Probs){

  //Set error for rounding issues
  double fuzz = 4.0*std::numeric_limits<double>::epsilon();

  //length of input vector
  int vecLen = Trials.ncol();

  //Setup things needed inside loop
  IntegerVector holdRow(vecLen), test(vecLen+2);
  NumericVector nppm(Probs.length()), h(Probs.length()), j(Probs.length());
  NumericMatrix RETURN(Trials.nrow(), Probs.length());


  //Loop over rows, calculate quantiles in each row
  for(int i=0; i<Trials.nrow(); i++){

    //subset matrix. I use several times
    holdRow = Trials.row(i);

    //sort vector
    std::sort(holdRow.begin(), holdRow.end());

    //Set vector with row, and an extra min/max. This is part of the formula,
    // It has been adapted here to avoid extra copies from push_front/bach
    test[0]=holdRow[0];
    test[seq(1,vecLen)] = holdRow;
    test[vecLen+1]=holdRow[vecLen-1];


    //1/3 and 2/3 are to match algorithm 8 in the rstats quantile function
    //Not really sure what else is happening. I took this from R quantile function
    nppm = 1.0/3.0 + Probs*(vecLen + 1.0/3.0);
    h = nppm + fuzz;
    j = floor(h);
    h = nppm - j;

    //Safety to set small values of h to 0. This does not handle negative values
    // like the r quantile does.
    h[h<fuzz] = 0;


    //calculate quantiles. Also not sure why or how this works.
    //No return value, set in place
    RETURN(i,_) = (1-h)*as<NumericVector>(test[j])+h*as<NumericVector>(test[j+1]);

  }//end loop over rows
  return(RETURN);
}//end function















