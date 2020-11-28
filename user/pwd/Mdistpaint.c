/* Geologic distance painting by plane-wave construction. */
/* Yunzhi Shi, March 2017 */
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

#include <math.h>
#include <rsf.h>
#include "predict.h"

void normalize(int n1, int n2, float **scan, float scaleFac);

int main (int argc, char *argv[])
{
    bool verb;
    int n1,n2,n3, i1,i2,i3, i0, order;
    float d2, eps, scaleFac, **u, **p, **fault, *trace;
    sf_file out, dip, flt=NULL, seed=NULL;

    sf_init(argc,argv);
    dip = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histfloat(dip,"d2",&d2)) d2=1.;
    n3 = sf_leftsize(dip,2);

    if (NULL != sf_getstring("seed")) {
	    seed = sf_input("seed");
    } else {
        sf_error("Required file [seed] not found");
    }
    if (NULL != sf_getstring("flt")) {
	    flt = sf_input("flt");
    } else {
        sf_error("Required file [flt] not found");
    }

    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("i0",&i0)) i0=0;
    /* reference trace */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    if (!sf_getfloat("faultscale",&scaleFac)) scaleFac=100.;
    /* Fault attribute scaling factor (0.0 ~ factor) */

    predict_init (n1, n2, eps*eps, order, 1, false);

    u = sf_floatalloc2(n1,n2);
    p = sf_floatalloc2(n1,n2);
    fault = sf_floatalloc2(n1,n2);
    trace = sf_floatalloc(n1);

    for (i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("cmp %d of %d;",i3+1,n3);
	sf_floatread(p[0],n1*n2,dip);

	sf_floatread(fault[0],n1*n2,flt);
    normalize(n1,n2,fault,scaleFac);

	sf_floatread(trace,n1,seed);

    // Initialize the reference trace
	for (i1=0; i1 < n1; i1++) u[i0][i1] = trace[i1];

    // Backward: accumulative compute geologic distance
	for (i2=i0-1; i2 >= 0; i2--) {
	    predict_step(false,false,trace,p[i2]);
	    for (i1=0; i1 < n1; i1++) u[i2][i1] = trace[i1];
        for (i1=0; i1<n1; i1++) trace[i1] += d2 * sqrt(1+p[i2][i1]*p[i2][i1])
                                           + fault[i2][i1];
	}

    // Forward: accumulative compute geologic distance
	for (i1=0; i1 < n1; i1++) trace[i1] = u[i0][i1];
	for (i2=i0+1; i2 < n2; i2++) {
	    predict_step(false,true,trace,p[i2-1]);
	    for (i1=0; i1 < n1; i1++) u[i2][i1] = trace[i1];
        for (i1=0; i1<n1; i1++) trace[i1] += d2 * sqrt(1+p[i2][i1]*p[i2][i1])
                                           + fault[i2][i1];
	}

	sf_floatwrite(u[0],n1*n2,out);
    }

    if (verb) sf_warning(".");

    exit (0);
}

void normalize(int n1, int n2, float **scan, float scaleFac) {
    int i1,i2;
    float min=FLT_MAX,max=0.0;

    for (i1=0; i1 < n1; i1++) {
        for (i2=0; i2 < n2; i2++) {
            if (scan[i2][i1] < min) min = scan[i2][i1];
            if (scan[i2][i1] > max) max = scan[i2][i1];
        }
    }
    if ((max-min) < 1e-6) sf_warning("WARNING: Input fault range < 1e-6.");

    for (i1=0; i1 < n1; i1++)
        for (i2=0; i2 < n2; i2++)
            scan[i2][i1] = (scan[i2][i1]-min) / (max-min) * scaleFac;
}
