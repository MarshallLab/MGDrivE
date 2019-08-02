#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************
 * Helpers
 *****************************************************************************/
/**************************************
 * Helper Constants
 *************************************/

/* limits for truncated distributions */
const static double inf_pos = std::numeric_limits<double>::infinity();
//const static double inf_neg = -std::numeric_limits<double>::infinity();

/**************************************
 * Helper Functions
 *************************************/

/* approximate equality */
template<typename T>
static bool approxEqual(T f1, T f2) {
  return (std::fabs(f1 - f2) <= std::numeric_limits<T>::epsilon() * fmax(std::fabs(f1), std::fabs(f2)));
}

/* truncated exponential distribution */
inline double dtruncExp(double x, double r, double a, double b){
  if(a >= b){
    Rcpp::stop("argument a is greater than or equal to b\n");
  }
  double scale = 1.0/r;
  double Ga = R::pexp(a,scale,true,false);
  double Gb = R::pexp(b,scale,true,false);
  if(approxEqual(Ga,Gb)){
    Rcpp::stop("Truncation interval is not inside the domain of the density function\n");
  }
  double density = R::dexp(x,scale,false) / (R::pexp(b,scale,true,false) - R::pexp(a,scale,true,false));
  return density;
}

/******************************************************************************
 * Kernels
 *****************************************************************************/
/**************************************
 * lognormal
 *************************************/

//' Calculate Lognormal Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}},
//' calculate a stochastic matrix where one step movement probabilities follow a lognormal density.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param meanlog log mean of \code{\link[stats]{Lognormal}} distribution
//' @param sdlog log standard deviation of \code{\link[stats]{Lognormal}} distribution
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate lognormal distribution over distances
//' #  mean and standard deviation are just for example
//' kernMat = calcLognormalKernel(distMat = distMat, meanlog = 100, sdlog = 10)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcLognormalKernel(const Rcpp::NumericMatrix& distMat,
                                        const double& meanlog, const double& sdlog){

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      kernMat(i,j) = R::dlnorm(distMat(i,j),meanlog,sdlog,false);
    }
    kernMat(i,_) = kernMat(i,_) / Rcpp::sum(kernMat(i,_)); /* normalize density */
  }

  return kernMat;
};

/**************************************
 * gamma
 *************************************/

//' Calculate Gamma Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}}, calculate a
//' stochastic matrix where one step movement probabilities follow a gamma density.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param shape shape parameter of \code{\link[stats]{GammaDist}} distribution
//' @param rate rate parameter of \code{\link[stats]{GammaDist}} distribution
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate gamma distribution over distances
//' #  shape and rate are just for example
//' kernMat = calcGammaKernel(distMat = distMat, shape = 1, rate = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcGammaKernel(const Rcpp::NumericMatrix& distMat, const double& shape, const double& rate){

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      kernMat(i,j) = R::dgamma(distMat(i,j),shape,rate,false);
    }
    kernMat(i,_) = kernMat(i,_) / Rcpp::sum(kernMat(i,_)); /* normalize density */
  }

  return kernMat;
};

/**************************************
 * exponential
 *************************************/

//' Calculate Exponential Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}}, calculate a
//' stochastic matrix where one step movement probabilities follow an exponential density.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param rate rate parameter of \code{\link[stats]{Exponential}} distribution
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate exponential distribution over distances
//' #  rate is just for example
//' kernMat = calcExpKernel(distMat = distMat, rate = 10)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcExpKernel(const Rcpp::NumericMatrix& distMat, const double& rate){

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);
  double scale = 1.0/rate;

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      kernMat(i,j) = R::dexp(distMat(i,j),scale,false);
    }
    kernMat(i,_) = kernMat(i,_) / Rcpp::sum(kernMat(i,_)); /* normalize density */
  }

  return kernMat;
};

/**************************************
 * hurdle exponential (point mass at zero + zero-truncated exponential distribution)
 *************************************/

//' Calculate Hurdle Exponential Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}}, calculate a
//' stochastic matrix where one step movement probabilities follow an zero-truncated
//' exponential density with a point mass at zero.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param rate rate parameter of \code{\link[stats]{Exponential}} distribution
//' @param pi point mass at zero
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate hurdle exponential distribution over distances
//' #  rate and point mass are just for example
//' kernMat = calcHurdleExpKernel(distMat = distMat, rate = 10, pi = 1000)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcHurdleExpKernel(const Rcpp::NumericMatrix& distMat, double rate, double pi){
  const double a = 1.0e-10; /* lower truncation bound */

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      if(i==j){
        kernMat(i,j) = 0;
      } else {
        kernMat(i,j) = dtruncExp(distMat(i,j),rate,a,inf_pos); /* truncated density */
      }
    }
    kernMat(i,_) = (kernMat(i,_) / Rcpp::sum(kernMat(i,_))) *(1-pi); /* normalize density */
    kernMat(i,i) = pi; /* point mass at zero */
  }

  return kernMat;
}
