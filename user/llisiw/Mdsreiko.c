/* Double square-root eikonal solver (2D) */
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

int main(int argc, char* argv[])
{
    bool velocity;
    int dim, i, n[SF_MAX_DIM], is, ns;
    int *f;
    float *th, *al;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *s, *t;
    float tau1, tau2, angle, thres;
    char key[6];
    sf_file in, out, flag, theta, alpha;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    /* read input dimension */
    dim = sf_filedims(in,n);

    ns = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input.",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	ns *= n[i];
    }
    if (dim < 3) {
	/* extend the third dimension for output (copy second dimension) */
	n[2] = n[1]; d[2] = d[1]; o[2] = o[1];
    }

    /* read input */
    s = sf_floatalloc(ns);
    sf_floatread(s,ns,in);

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */

    /* convert to slowness squared */
    if (velocity) {
	for (is=0; is < ns; is++)
	    s[is] = 1./s[is]*1./s[is];
    }

    if (!sf_getfloat("tau1",&tau1)) tau1=1.e-3;
    /* tau1 */

    if (!sf_getfloat("tau2",&tau2)) tau2=1.;
    /* tau2 */

    if (!sf_getfloat("angle",&angle)) angle=5.;
    /* angle (degree) */

    angle = tan(angle/180.*3.1416);

    if (!sf_getfloat("thres",&thres)) thres=0.1;
    /* threshold (percentage) */

    if (NULL != sf_getstring("flag")) {
	flag = sf_output("flag");
	sf_settype(flag,SF_INT);
	sf_putint(flag,"n3",n[1]);
	sf_putfloat(flag,"d3",d[1]);
	sf_putfloat(flag,"o3",o[1]);
	f = sf_intalloc(ns*n[1]);
    } else {
	flag = NULL;
	f = NULL;
    }

    if (NULL != sf_getstring("theta")) {
	theta = sf_output("theta");
	sf_putint(theta,"n3",n[1]);
	sf_putfloat(theta,"d3",d[1]);
	sf_putfloat(theta,"o3",o[1]);
	th = sf_floatalloc(ns*n[1]);
    } else {
	theta = NULL;
	th = NULL;
    }

    if (NULL != sf_getstring("alpha")) {
	alpha = sf_output("alpha");
	sf_putint(alpha,"n3",n[1]);
	sf_putfloat(alpha,"d3",d[1]);
	sf_putfloat(alpha,"o3",o[1]);
	al = sf_floatalloc(ns*n[1]);
    } else {
	alpha = NULL;
	al = NULL;
    }

    /* allocate memory for output */
    t = sf_floatalloc(ns*n[1]);

    /* initialize */
    dsreiko_init(n,o,d,tau1,tau2,angle,thres);

    /* compute */
    dsreiko_fastmarch(t,s,f,th,al);

    /* mirror */
    dsreiko_mirror(t);

    if (flag != NULL) sf_intwrite(f,ns*n[1],flag);
    if (theta != NULL) sf_floatwrite(th,ns*n[1],theta);
    if (alpha != NULL) sf_floatwrite(al,ns*n[1],alpha);

    /* write output dimension */
    sf_putint(out,"n3",n[1]);
    sf_putfloat(out,"d3",d[1]);
    sf_putfloat(out,"o3",o[1]);

    /* write output */
    sf_floatwrite(t,ns*n[1],out);

    exit(0);
}
