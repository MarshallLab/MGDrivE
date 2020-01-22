#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************
 * Helper Functions
 *****************************************************************************/

/* Convert degrees to radians */
inline double toRad(const double& deg) {
  return deg * 0.0174532925199433;
}

/* converts radians to degrees */
inline double toDeg(const double& rad) {
  return rad * 57.2957795130823 ;
}

/******************************************************************************
 * Distance Functions
 *****************************************************************************/

/* These functions come from the geosphere package, version 1.5-7
 * They are pulled from src/dist.c file, and converted to c++
 *
 * The original functions (by Sean) will be commented out above the ones used
 * for the package. The originals are less clear to read, but generally the same.
 */

/******************************************************************************
 * Cosine Distance
 *****************************************************************************/

// double distCos(double lon1, double lat1, double lon2, double lat2, double r) {
// 	double cosd;
// 	lon1 = toRad(lon1);
// 	lon2 = toRad(lon2);
// 	lat1 = toRad(lat1);
// 	lat2 = toRad(lat2);
// 	cosd = std::sin(lat1) * std::sin(lat2) + std::cos(lat1) * std::cos(lat2) * std::cos(lon1-lon2);
// 	return std::acos(cosd) * r;
// }


//' Calculate Geodesic Distance - Cosine Method
//'
//' This function calculates geodesic distance using the cosine method.
//'
//' @param latLongs Two column matrix of latitudes/longitudes
//' @param r Earth radius. Default is WGS-84 radius
//'
//' @examples
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # cosine distance formula
//' distMat = calcCos(latLongs = latLong)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcCos(const Rcpp::NumericMatrix& latLongs, const double& r = 6378137){

  // latLongs = 2 column matrix of latitudes/longitudes in degrees
  // r = equatorial radius of earth in meters from WGS-84

  ////////////////////
  // Setup objects
  ////////////////////
  size_t n(latLongs.nrow());
  Rcpp::NumericMatrix retMat(n, n);
  double lon1, lat1, lon2, lat2, cosd;


  ////////////////////
  // Loop over all points, calculating distance
  ////////////////////
  for(size_t row = 1; row < n; ++row){

    // set lat/long for first point
    lon1 = toRad(latLongs(row,1));
    lat1 = toRad(latLongs(row,0));


    for(size_t col = 0; col < row; ++col){

      // set lat/long for second point
      lon2 = toRad(latLongs(col,1));
      lat2 = toRad(latLongs(col,0));

      cosd = std::sin(lat1) * std::sin(lat2) + std::cos(lat1) *
                        std::cos(lat2) * std::cos(lon1-lon2);

      // fill return matrix
      retMat(row,col) = retMat(col,row) = std::acos(cosd) * r;

    }// end loop over columns
  }// end loop over rows

  // return matrix
  return(retMat);

}// end function

/******************************************************************************
 * Haversine Distance
 *****************************************************************************/

/********************
 * Seans original Haversine function, with changed radius to reflect WGS-84
 ******************/
// double gcd_hf(const double& long1, const double& lat1, const double& long2, const double& lat2){
//   double deltaLong = (long2 - long1);
//   double deltaLat = (lat2 - lat1);
//   double a = std::pow(std::sin(deltaLat/2),2) + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(deltaLong/2),2);
//   double sqrtA = std::sqrt(a);
//   double c = 2 * std::asin(fmin(1.0,sqrtA));
//   double d = c * 6378137;
//   //double d = (6371.0 * c)*1000.0;
//   return d;
// };
//
// Rcpp::NumericMatrix calc_Haversine(const Rcpp::NumericMatrix& latlongs){
//   size_t n = latlongs.nrow();
//   Rcpp::NumericMatrix zz = Rcpp::NumericMatrix(n,n);
//   for(size_t i=0; i<n; i++){
//     for(size_t j=0; j<n; j++){
//       zz(i,j) = gcd_hf(toRad(latlongs(i,1)), toRad(latlongs(i,0)),toRad(latlongs(j,1)), toRad(latlongs(j,0)));
//     }
//   }
//   return zz;
// };
//
// dist.c function
// double distHav(double lon1, double lat1, double lon2, double lat2, double r) {
//
// 	double dLat, dLon, a;
//
// 	lon1 = toRad(lon1);
// 	lon2 = toRad(lon2);
// 	lat1 = toRad(lat1);
// 	lat2 = toRad(lat2);
//
// 	dLat = lat2-lat1;
// 	dLon = lon2-lon1;
// 	a = sin(dLat/2.) * sin(dLat/2.) + cos(lat1) * cos(lat2) * sin(dLon/2.) * sin(dLon/2.);
// 	return 2. * atan2(sqrt(a), sqrt(1.-a)) * r;
// }


//' Calculate Geodesic Distance - Haversine Method
//'
//' This function calculates geodesic distance using the Haversine method.
//'
//' @param latLongs Two column matrix of latitudes/longitudes
//' @param r Earth radius. Default is WGS-84 radius
//'
//' @examples
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Haversine distance formula
//' distMat = calcHaversine(latLongs = latLong)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcHaversine(const Rcpp::NumericMatrix& latLongs, const double& r = 6378137){

  // latLongs = 2 column matrix of latitudes/longitudes in degrees
  // r = equatorial radius of earth in meters from WGS-84

  ////////////////////
  // Setup objects
  ////////////////////
  size_t n(latLongs.nrow());
  Rcpp::NumericMatrix retMat(n, n);
  double lon1, lat1, lon2, lat2, dLat, dLon, a;


  ////////////////////
  // Loop over all points, calculating distance
  ////////////////////
  for(size_t row = 1; row < n; ++row){

    // set lat/long for first point
    lon1 = toRad(latLongs(row,1));
    lat1 = toRad(latLongs(row,0));


    for(size_t col = 0; col < row; ++col){

      // set lat/long for second point
      lon2 = toRad(latLongs(col,1));
      lat2 = toRad(latLongs(col,0));


      dLat = lat2-lat1;
	    dLon = lon2-lon1;
	    a = std::sin(dLat/2.0) * std::sin(dLat/2.0)
	    	+ std::cos(lat1) * std::cos(lat2) * std::sin(dLon/2.0) * std::sin(dLon/2.0);


      // fill return matrix
      retMat(row,col) = retMat(col,row) = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0-a)) * r;

    }// end loop over columns
  }// end loop over rows

  // return matrix
  return(retMat);

}// end function

/******************************************************************************
 * Vincenty Sphere Distance
 *****************************************************************************/

// double distVinSph(double lon1, double lat1, double lon2, double lat2, double r) {
//
// 	double x, x1, x2, y;
//
// 	lon1 = toRad(lon1);
// 	lon2 = toRad(lon2);
// 	lat1 = toRad(lat1);
// 	lat2 = toRad(lat2);
//
// 	// These are incorrect.
// 	// x1 = sqrt(cos(lat2) * sin(lon1-lon2));
// 	// x2 = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2);
// 	// x = x1*x1 + x2*x2;
//
// 	// replacement x things
// 	x1 = cos(lat2) * sin(lon1-lon2);
//   x2 = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1-lon2);
//   x = sqrt(x1*x1 + x2*x2);
// 	y = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2);
//
// 	return r * atan2(x, y);
// }

//' Calculate Geodesic Distance - Vincenty Sphere Method
//'
//' This function calculates geodesic distance using the Vincenty sphere method.
//'
//' @param latLongs Two column matrix of latitudes/longitudes
//' @param r Earth radius. Default is WGS-84 radius
//'
//' @examples
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Sphere  distance formula
//' distMat = calcVinSph(latLongs = latLong)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcVinSph(const Rcpp::NumericMatrix& latLongs, const double& r = 6378137){

  // latLongs = 2 column matrix of latitudes/longitudes in degrees
  // r = equatorial radius of earth in meters from WGS-84

  ////////////////////
  // Setup objects
  ////////////////////
  size_t n(latLongs.nrow());
  Rcpp::NumericMatrix retMat(n, n);
  double lon1, lat1, lon2, lat2, x, x1, x2, y;


  ////////////////////
  // Loop over all points, calculating distance
  ////////////////////
  for(size_t row = 1; row < n; ++row){

    // set lat/long for first point
    lon1 = toRad(latLongs(row,1));
    lat1 = toRad(latLongs(row,0));

    for(size_t col = 0; col < row; ++col){

      // set lat/long for second point
      lon2 = toRad(latLongs(col,1));
      lat2 = toRad(latLongs(col,0));

      x1 = std::cos(lat2) * std::sin(lon1-lon2);
      x2 = std::cos(lat1) * std::sin(lat2) - std::sin(lat1) * std::cos(lat2) * std::cos(lon1-lon2);
      x = std::sqrt(x1*x1 + x2*x2);
    	y = std::sin(lat1) * std::sin(lat2) + std::cos(lat1) * std::cos(lat2) * std::cos(lon1-lon2);

      // fill return matrix
      retMat(row,col) = retMat(col,row) = r * atan2(x, y);

    }// end loop over columns
  }// end loop over rows

  // return matrix
  return(retMat);

}// end function

/******************************************************************************
 * Vincenty Ellipse Distance
 *****************************************************************************/

// double distVinEll(double lon1, double lat1, double lon2, double lat2, double a, double b, double f) {
//   /*  Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2009
//    Calculate geodesic distance (in m) between two points specified by latitude/longitude
//    (in numeric degrees) using Vincenty inverse formula for ellipsoids
//    based on source http://www.movable-type.co.uk/scripts/latlong-vincenty.html (c) 2002-2009 Chris Veness
//    */
//   double L, U1, U2, sinU1, cosU1, sinU2, cosU2, lambda, sinLambda, cosLambda, sinSigma,
//   cosSigma, sigma, sinAlpha, cosSqAlpha, cos2SigmaM, C, lambdaP, uSq, A, B, deltaSigma;
//
//   int iterLimit, cont;
//
//   if ((lon1 == lon2) & (lat1 == lat2))  {
//
//     return 0.0;
//
//   } else {
//
//     lon1 = toRad(lon1);
//     lon2 = toRad(lon2);
//     lat1 = toRad(lat1);
//     lat2 = toRad(lat2);
//
//     L = (lon2-lon1);
//     U1 = std::atan((1.-f) * std::tan(lat1));
//     U2 = std::atan((1.-f) * tan(lat2));
//     sinU1 = std::sin(U1);
//     cosU1 = std::cos(U1);
//     sinU2 = std::sin(U2);
//     cosU2 = std::cos(U2);
//     lambda = L;
//     iterLimit = 100;
//     cont = 1;
//
//     while (cont) {
//       sinLambda = std::sin(lambda);
//       cosLambda = std::cos(lambda);
//       sinSigma = std::sqrt((cosU2*sinLambda) * (cosU2*sinLambda) + (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda));
//
//       cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
//       sigma = std::atan2(sinSigma, cosSigma);
//       sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
//       cosSqAlpha = 1. - sinAlpha*sinAlpha;
//       cos2SigmaM = cosSigma - 2.*sinU1*sinU2/cosSqAlpha;
//
//       // if (Rcpp::is_na(cos2SigmaM)) {
//       //   cos2SigmaM = 0.;  // equatorial line: cosSqAlpha=0 (?6)
//       // }
//
//       if (Rcpp::traits::is_na<REALSXP>(cos2SigmaM)) {
//         cos2SigmaM = 0.;  // equatorial line: cosSqAlpha=0 (?6)
//       }
//
//       C = f/16.*cosSqAlpha*(4.+f*(4.-3.*cosSqAlpha));
//       lambdaP = lambda;
//       lambda = L + (1.-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1.+2.*cos2SigmaM*cos2SigmaM)));
//
//       iterLimit = iterLimit - 1;
//       cont = (std::abs(lambda-lambdaP) > 1e-12 && iterLimit > 0);
//     }
//
//     if (iterLimit==0) {
//       return R_NaReal;  // formula failed to converge
//     } else {
//       uSq = cosSqAlpha * (a*a - b*b) / (b*b);
//       A = 1.0 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175.0*uSq)));
//       B = uSq/1024.0 * (256.0 + uSq*(-128.0 + uSq*(74.0-47.0*uSq)));
//       deltaSigma = B*sinSigma*(cos2SigmaM+B/4.*(cosSigma*(-1.+2.*cos2SigmaM*cos2SigmaM)- B/6.*cos2SigmaM*(-3.+4.*sinSigma*sinSigma)*(-3.+4.*cos2SigmaM*cos2SigmaM)));
//       return  b*A*(sigma-deltaSigma);
//     }
//   }
// }


//' Calculate Geodesic Distance - Vincenty Ellipsoid Method
//'
//' This function calculates geodesic distance using the original Vincenty Ellipsoid method.
//'
//' @param latLongs Two column matrix of latitudes/longitudes
//' @param a Equatorial radius of the earth, default is WGS-84 radius
//' @param b Polar radius of the earth, default is WGS-84 radius
//' @param f Flattening or inverse eccentricity, default eccentricity is WGS-84
//' @param eps Convergence criteria
//' @param iter Maximum number of iterations to attempt convergence
//'
//' @examples
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcVinEll(const Rcpp::NumericMatrix& latLongs, const double& a = 6378137,
                                  const double& b = 6356752.3142, const double& f = 1.0/298.257223563,
                                  const double& eps = 1e-12, const double& iter = 100) {

  // latLongs = 2 column matrix of latitudes/longitudes in degrees
  // a = equatorial radius of earth in meters from WGS-84
  // b = polar radius of earth in meters from WGS-84
  // f = flattening, 1/ eccentricity of the earth from WGS-84
  // eps = convergence criteria
  // iter = how many iterations to test and hope it converges?

  /*  Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2009
   Calculate geodesic distance (in m) between two points specified by latitude/longitude
   (in numeric degrees) using Vincenty inverse formula for ellipsoids
   based on source http://www.movable-type.co.uk/scripts/latlong-vincenty.html (c) 2002-2009 Chris Veness
   */


  ////////////////////
  // Setup objects
  ////////////////////
  size_t n(latLongs.nrow());
  Rcpp::NumericMatrix retMat(n, n);

  double L, U1, U2, sinU1, cosU1, sinU2, cosU2, lambda, sinLambda, cosLambda, sinSigma,
  cosSigma, sigma, sinAlpha, cosSqAlpha, cos2SigmaM, C, lambdaP, uSq, A, B, deltaSigma,
  lon1, lat1, lon2, lat2;

  int iterLimit;
  bool cont;


  ////////////////////
  // Loop over all points, calculating distance
  ////////////////////
  for(size_t row = 1; row < n; ++row){

    // set lat/long for first point
    lon1 = toRad(latLongs(row,1));
    lat1 = toRad(latLongs(row,0));


    for(size_t col = 0; col < row; ++col){

      // set lat/long for second point
      lon2 = toRad(latLongs(col,1));
      lat2 = toRad(latLongs(col,0));


      // calulate initial things or things that are reused
      L = (lon2-lon1);
      U1 = std::atan((1.0-f) * std::tan(lat1));
      U2 = std::atan((1.0-f) * tan(lat2));
      sinU1 = std::sin(U1);
      cosU1 = std::cos(U1);
      sinU2 = std::sin(U2);
      cosU2 = std::cos(U2);
      lambda = L;
      iterLimit = iter;
      cont = true;

      // convergence algorithm! (algorithm? kinda shitty for that...)
      while (cont) {
        sinLambda = std::sin(lambda);
        cosLambda = std::cos(lambda);
        sinSigma = std::sqrt((cosU2*sinLambda) * (cosU2*sinLambda) +
                    (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda));

        cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
        sigma = std::atan2(sinSigma, cosSigma);
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1.0 - sinAlpha*sinAlpha;
        cos2SigmaM = cosSigma - 2.0*sinU1*sinU2/cosSqAlpha;

        if (Rcpp::traits::is_na<REALSXP>(cos2SigmaM)) {
          cos2SigmaM = 0.;  // equatorial line: cosSqAlpha=0 (?6)
        }

        C = f/16.0*cosSqAlpha * (4.0 + f*(4.0 - 3.0*cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1.0-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1.0 +2.0*cos2SigmaM*cos2SigmaM)));

        // check if converged or out of iterations
        iterLimit--;
        cont = (std::abs(lambda-lambdaP) > eps && iterLimit > 0);
      }// end while loop

      // check if converged or not, and fill return matrix
      if (iterLimit==0) {
        retMat(row,col) = retMat(col,row) =  R_NaReal;  // formula failed to converge
      } else {
        uSq = cosSqAlpha * (a*a - b*b) / (b*b);
        A = 1.0 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175.0*uSq)));
        B = uSq/1024.0 * (256.0 + uSq*(-128.0 + uSq*(74.0-47.0*uSq)));
        deltaSigma = B*sinSigma*(cos2SigmaM + B/4.0 * (cosSigma*(-1.0 + 2.0*cos2SigmaM*cos2SigmaM)
                                                         - B/6.0 * cos2SigmaM*(-3.0 + 4.0*sinSigma*sinSigma)
                                                         *(-3.0 + 4.0*cos2SigmaM*cos2SigmaM)));
        // fill return matrix
        retMat(row,col) = retMat(col,row) = b*A*(sigma-deltaSigma);
      }// end check and fill

    }// end loop over columns
  }// end loop over rows

  // return matrix
  return(retMat);

}// end function


