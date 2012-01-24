/* Double square-root eikonal solver interface */
/*
  Copyright (C) 2011 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

#include "dsreiko.h"

const double tol = 1.e-6;

void ferrari(double a, double b, double c, double d, double e /* coefficients */,
	     double *root /* root */)
/*< quartic solve (Ferrari's method) >*/
{
    /* a*x^4 + b*x^3 + c*x^2 + d*x + e = 0. */
    /* http://en.wikipedia.org/wiki/Quartic_function */

    double alpha, beta, gama;
    double y, P, Q, R, U, W;

    /* first step */
    alpha = -3./8.*pow(b,2.)/pow(a,2.)+c/a;
    beta  = 1./8.*pow(b,3.)/pow(a,3.)-1./2.*b*c/pow(a,2.)+d/a;
    gama  = -3./256.*pow(b,4.)/pow(a,4.)+1./16.*c*pow(b,2.)/pow(a,3.)-1./4.*b*d/pow(a,2.)+e/a;

    /* second step */
    if (fabs(beta) <= tol) {
	root[0] = sf_cmplx(-1./4.*b/a+sqrt((-alpha+sqrt(pow(alpha,2.)-4.*gama))/2.),0.);
	root[1] = sf_cmplx(-1./4.*b/a-sqrt((-alpha+sqrt(pow(alpha,2.)-4.*gama))/2.),0.);
	root[2] = sf_cmplx(-1./4.*b/a+sqrt((-alpha-sqrt(pow(alpha,2.)-4.*gama))/2.),0.);
	root[3] = sf_cmplx(-1./4.*b/a-sqrt((-alpha-sqrt(pow(alpha,2.)-4.*gama))/2.),0.);

	return();
    }

    /* third step */
    P = -1./12.*pow(alpha,2.)-gama;
    Q = -1./108.*pow(alpha,3.)+1./3.*alpha*gama-1./8.*pow(beta,2.);
    R = -1./2.*Q+sqrt(1./4.*pow(Q,2.)+1./27.*pow(P,3.));
    /* R = -1./2.*Q-sqrt(1./4.*pow(Q,2.)+1./27.*pow(P,3.)); */

    U = pow(R,1./3.);

    if (fabs(U) <= tol)
	y = -5./6.*alpha+U-pow(Q,1./3.);
    else
	y = -5./6.*alpha+U-1./3.*P/U;

    W = sqrt(alpha+2.*y);

    /* final step */
    root[0] = -1./4.*b/a+(W+sqrt(-(3.*alpha+2.*y+2.*beta/W)))/2.;
    root[1] = -1./4.*b/a+(-W+sqrt(-(3.*alpha+2.*y-2.*beta/W)))/2.;
    root[2] = -1./4.*b/a+(W-sqrt(-(3.*alpha+2.*y+2.*beta/W)))/2.;
    root[3] = -1./4.*b/a+(-W-sqrt(-(3.*alpha+2.*y-2.*beta/W)))/2.;
}
