#ifndef __SEAM_GAUSS_WAVELET__
#define __SEAM_GAUSS_WAVELET__

#include "cstd.h"
#include "utils.h"

/** returns array containing Ricker wavelet, symmetric about midpoint,
    computes length 2*iw+1 to truncate sampling at approximately
    single precision round-off error (i.e. O(10^7)). Scaled to unit
    amplitude at t=0 (sample iw+1).
*/
float * getrick(int * iw, float dt, float fpeak);

/** returns array containing Gaussian wavelet, symmetric about
    midpoint, of length 2*iw+1 to truncate sampling at approximately
    single precision round-off error (i.e. O(10^7)). Scaled so that
    2nd derivative is Ricker with unit amplitude at t=0 (sample iw+1).
*/
float * getgauss(int * iw, float dt, float fpeak);
ireal * igetgauss(int * iw, ireal dt, ireal fpeak);
/** returns array containing derivative of Gaussian wavelet, symmetric about
    midpoint, of length 2*iw+1 to truncate sampling at approximately
    single precision round-off error (i.e. O(10^7)). Scaled so that
    1st derivative is Ricker with unit amplitude at t=0 (sample iw+1).
*/
float * getdgauss(int * iw, float dt, float fpeak);
ireal * igetdgauss(int * iw, ireal dt, ireal fpeak);


/** returns value of gaussian at a point, scaled by -(2 pi^2 fpeak^2)^-1 */
float compgauss(float r, float fpeak);

/** returns value of gaussian derivative at a point scaled by -(2 pi^2 fpeak^2)^-1  */
float compdgauss(float r, float fpeak);

/** returns value of gaussian 2ndderivative at a point scaled by -(2 pi^2 fpeak^2)^-1 - same as normalized Ricker */
  float comprick(float r, float fpeak);

#endif
