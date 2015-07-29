/* Fast marching for image rays */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include "irays.h"

int main(int argc, char* argv[])
{
    bool velocity;
    int dim, i, n[3], it, nt, order, j;
    float d[3], o[3], thres;
    float *s, *t0, *x0;
    int *f0;
    char key[3];
    sf_file in, out, ot0=NULL, ox0=NULL, of0=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* read input dimension */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i];
    }
    if (dim < 3) {
	n[2] = 1; d[2] = d[1]; o[2] = o[1];
    }

    /* read initial guess */
    s = sf_floatalloc(nt);
    sf_floatread(s,nt,in);

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* y, inputs are velocity / n, slowness-squared */

    if (velocity) {
	for (it=0; it < nt; it++) {
	    s[it] = 1./s[it]*1./s[it];
	}
    }

    /* allocate temporary memory */
    t0 = sf_floatalloc(nt);
    x0 = sf_floatalloc(nt);
    f0 = sf_intalloc(nt);

    /* output transformation matrix */
    if (NULL != sf_getstring("t0"))
	ot0 = sf_output("t0");
    if (NULL != sf_getstring("x0"))
	ox0 = sf_output("x0");

    /* output upwind neighbor */
    if (NULL != sf_getstring("f0")) {
	of0 = sf_output("f0");
	sf_settype(of0,SF_INT);
    }

    if (!sf_getint("order",&order)) order=1;
    /* fastmarching accuracy order */

    if (!sf_getfloat("thres",&thres)) thres=10.;
    /* thresholding for caustics */

    /* initialization */
    fastmarch_init(n,o,d,order);

    /* fastmarch */
    fastmarch(t0,x0,f0,s);

    /* write output */
    sf_floatwrite(s,nt,out);
    if (NULL!=ot0) sf_floatwrite(t0,nt,ot0);
    if (NULL!=ox0) sf_floatwrite(x0,nt,ox0);

    /* caustic region (2D) */
    for (i=0; i < n[0]; i++) {
	for (j=0; j < n[1]; j++) {
	    if (j > 0) {
		if (x0[j*n[0]+i] <= x0[(j-1)*n[0]+i]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
		if ((x0[j*n[0]+i]-x0[(j-1)*n[0]+i]) > thres*d[1]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
	    }
	    if (j < n[1]-1) {
		if (x0[(j+1)*n[0]+i] <= x0[j*n[0]+i]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
		if ((x0[(j+1)*n[0]+i]-x0[j*n[0]+i]) > thres*d[1]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
	    }
	}
    }

    /* write flag */
    if (NULL!=of0) sf_intwrite(f0,nt,of0);

    exit(0);
}
