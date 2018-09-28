/* Roatation with Interpolation from a regular grid in 2-D. */
/*
  Copyright (C) 2004 University of Texas at Austin

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
#include <math.h>
#include "shprefilter.h"

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2, k1, k2;
    float x1, x2, c1, c2, cosa, sina, angle;
    float **orig, **rotd;
    char *interp;
    float d1,d2;
    sf_file inp, out;

    /* initialize */
    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    /* get dimensions from input */
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in inp");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in inp");
    if (!sf_histfloat(inp,"d1",&d1)) sf_error("No d1= in inp");
    if (!sf_histfloat(inp,"d2",&d2)) sf_error("No d2= in inp");

    /* get parameters from command line */
    if (!sf_getfloat("angle",&angle)) angle=90.;
    /* rotation angle */

    if (NULL == (interp = sf_getstring("interp"))) 
	interp="nearest";
    /* [n,l,c] interpolation type */

    /* convert degrees to radians */
    angle *= SF_PI/180.;
    cosa = cosf(angle);
    sina = sinf(angle);

    orig = sf_floatalloc2(n1,n2);
    rotd = sf_floatalloc2(n1,n2);

    /* read data */
    sf_floatread(orig[0],n1*n2,inp);

    /* central point */
    c1 = (n1-1)*0.5;
    c2 = (n2-1)*0.5;
    
/*     shifted linear interpolation prefilter */
    if (interp[0] == 's') {
    	shprefilter2d(n1,n2,orig); 
    }
    
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {

	    /* rotated coordinates */
	    x1 = c1+(i1-c1)*cosa-(i2-c2)*sina;
	    x2 = c2+(i1-c1)*sina+(i2-c2)*cosa;

	    /* nearest neighbor */
	    if(interp[0] != 's') {
	    	k1 = floorf(x1); x1 -= k1;
	    	k2 = floorf(x2); x2 -= k2;
	    } else { /*shift data (shifted interpolation)*/
	    	x1 -= 0.21*d1;
	    	x2 -= 0.21*d2;
	    	//sf_warning("td %f %f x %f %f ",tau*d1,tau*d2,x1,x2);
	    	k1 = floorf(x1); x1 -= k1;
	    	k2 = floorf(x2); x2 -= k2;
	    	//sf_warning("x %f %f k %f %f",x1,x2,k1,k2);
	    }

	    switch(interp[0]) {
		case 'n': /* nearest neighbor */
		    if (x1 > 0.5) k1++;
		    if (x2 > 0.5) k2++;
		    if (k1 >=0 && k1 < n1 &&
			k2 >=0 && k2 < n2) {
			rotd[i2][i1] = orig[k2][k1];
		    } else {
			rotd[i2][i1] = 0.;
		    }
		    break;
		case 'l': /* bilinear */
		    if (k1 >=0 && k1 < n1-1 &&
			k2 >=0 && k2 < n2-1) {
			rotd[i2][i1] = 
			    (1.-x1)*(1.-x2)*orig[k2][k1]   +
			    x1     *(1.-x2)*orig[k2][k1+1] +
			    (1.-x1)*x2     *orig[k2+1][k1] +
			    x1     *x2     *orig[k2+1][k1+1];
			    //sf_warning("%f",orig[k2][k1]);
		    } else {
			rotd[i2][i1] = 0.;
		    }
		    break;
		case 'c': /* cubic convolution */
		    if (k1 >=1 && k1 < n1-2 &&
			k2 >=1 && k2 < n2-2) {
			rotd[i2][i1] = 
(-0.5*x1*(1-x1)*(1-x1))*(-0.5*x2*(1-x2)*(1-x2))*orig[k2-1][k1-1]	+
(-0.5*x1*(1-x1)*(1-x1))*(0.5*(1-x2)*(2+2*x2-3*x2*x2))*orig[k2][k1-1]	+
(-0.5*x1*(1-x1)*(1-x1))*(0.5*x2*(1+4*x2-3*x2*x2))*orig[k2+1][k1-1]	+
(-0.5*x1*(1-x1)*(1-x1))*(-0.5*(1-x2)*x2*x2)*orig[k2+2][k1-1]		+
(0.5*(1-x1)*(2+2*x1-3*x1*x1))*(-0.5*x2*(1-x2)*(1-x2))*orig[k2-1][k1]	+
(0.5*(1-x1)*(2+2*x1-3*x1*x1))*(0.5*(1-x2)*(2+2*x2-3*x2*x2))*orig[k2][k1]+
(0.5*(1-x1)*(2+2*x1-3*x1*x1))*(0.5*x2*(1+4*x2-3*x2*x2))*orig[k2+1][k1]	+
(0.5*(1-x1)*(2+2*x1-3*x1*x1))*(-0.5*(1-x2)*x2*x2)*orig[k2+2][k1]	+
(0.5*x1*(1+4*x1-3*x1*x1))*(-0.5*x2*(1-x2)*(1-x2))*orig[k2-1][k1+1]	+
(0.5*x1*(1+4*x1-3*x1*x1))*(0.5*(1-x2)*(2+2*x2-3*x2*x2))*orig[k2][k1+1]	+
(0.5*x1*(1+4*x1-3*x1*x1))*(0.5*x2*(1+4*x2-3*x2*x2))*orig[k2+1][k1+1]	+
(0.5*x1*(1+4*x1-3*x1*x1))*(-0.5*(1-x2)*x2*x2)*orig[k2+2][k1+1]		+
(-0.5*(1-x1)*x1*x1)*(-0.5*x2*(1-x2)*(1-x2))*orig[k2-1][k1+2]		+
(-0.5*(1-x1)*x1*x1)*(0.5*(1-x2)*(2+2*x2-3*x2*x2))*orig[k2][k1+2]	+
(-0.5*(1-x1)*x1*x1)*(0.5*x2*(1+4*x2-3*x2*x2))*orig[k2+1][k1+2]		+
(-0.5*(1-x1)*x1*x1)*(-0.5*(1-x2)*x2*x2)*orig[k2+1][k1+2];
		    } else {
			rotd[i2][i1] = 0.;
		    }
		    break;
		case 's':  {/* shifted linear interpolation (bilinear + shifted + prefilter)*/
		    if (k1 >=0 && k1 < n1-1 &&
			k2 >=0 && k2 < n2-1) {
			rotd[i2][i1] = 
			    (1.-x1)*(1.-x2)*orig[k2][k1]   +
			    x1     *(1.-x2)*orig[k2][k1+1] +
			    (1.-x1)*x2     *orig[k2+1][k1] +
			    x1     *x2     *orig[k2+1][k1+1];
		    } else {
			rotd[i2][i1] = 0.;
		    }
		    break;
		}
		default:
		    sf_error("Unknown interpolation %s",
			     interp);
		    break;
	    }
	}
    }

    /* write result */
    sf_floatwrite(rotd[0],n1*n2,out);

    exit(0);
}
