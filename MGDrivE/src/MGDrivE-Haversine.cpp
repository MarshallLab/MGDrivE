#include <Rcpp.h>
#include <cmath>
#include <iostream>
using namespace Rcpp;

/* Convert degrees to radians */
inline double deg2rad(const double& deg){return deg*M_PI/180.0;};

/* Earth mean radius [km] */
const static double R_earth = 6371.0;

/* Calculates the geodesic distance between two points specified by
 * radian latitude/longitude using the Haversine formula
 * Ouputs distance between sites 1 and 2 as meters
*/

inline double gcd_hf(const double& long1, const double& lat1, const double& long2, const double& lat2){
  double deltaLong = (long2 - long1);
  double deltaLat = (lat2 - lat1);
  double a = std::pow(std::sin(deltaLat/2),2) + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(deltaLong/2),2);
  double sqrtA = std::sqrt(a);
  double c = 2 * std::asin(fmin(1.0,sqrtA));
  double d = (R_earth * c)*1000.0;
  return d;
};

/* Fxn to calculate matrix of distances between each two sites
 * INPUT: a matrix in which longs are in first column and lats in second column
 * OUTPUT: a distance matrix between all pairwise sites
 * Output distances are in meters
*/

//' Calculate Haversine Distance
//'
//' Calculate the great-circle distance (Haversine distance) between sets of longitude - latitude points, see \url{https://en.wikipedia.org/wiki/Haversine_formula}
//'
//' @param latlongs numeric matrix where first column is vector of latitudes and second column is vector of longitudes
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_haversine(const Rcpp::NumericMatrix& latlongs){
  size_t n = latlongs.nrow();
  Rcpp::NumericMatrix zz = Rcpp::NumericMatrix(n,n);
  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      zz(i,j) = gcd_hf(deg2rad(latlongs(i,1)), deg2rad(latlongs(i,0)),deg2rad(latlongs(j,1)), deg2rad(latlongs(j,0)));
    }
  }
  return zz;
};
