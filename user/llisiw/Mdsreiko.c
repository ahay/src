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
    bool velocity, causal;
    int dim, i, n[SF_MAX_DIM], nm;
    long nt;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *s, *t;
    int *f, nloop, *dp;
    float thres, *al, tol;
    char key[6];
    sf_file in, out, mask, flag, alpha;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    /* read input dimension */
    dim = sf_filedims(in,n);

    nm = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input.",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nm *= n[i];
    }
    if (dim < 3) {
	/* extend the third dimension for output (copy second dimension) */
	n[2] = n[1]; d[2] = d[1]; o[2] = o[1];
    }

    /* read input */
    s = sf_floatalloc(nm);
    sf_floatread(s,nm,in);

    /* write output dimension */
    sf_putint(out,"n3",n[2]);
    sf_putfloat(out,"d3",d[2]);
    sf_putfloat(out,"o3",o[2]);

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */

    /* convert to slowness squared */
    if (velocity) {
	for (i=0; i < nm; i++)
	    s[i] = 1./s[i]*1./s[i];
    }

    if (!sf_getfloat("thres",&thres)) thres=5.e-5;
    /* threshold (percentage) */

    if (!sf_getfloat("tol",&tol)) tol=1.e-3;
    /* tolerance for bisection root-search */

    if (!sf_getint("nloop",&nloop)) nloop=10;
    /* number of bisection root-search */

    nt = (long) nm*n[2];

    /* read receiver mask */    
    if (NULL == sf_getstring("mask")) {
	mask = NULL;
	dp = NULL;
    } else {
	mask = sf_input("mask");
	dp = sf_intalloc(n[1]*n[2]);
	sf_intread(dp,n[1]*n[2],mask);
	sf_fileclose(mask);
    }

    /* output DSR branch flag */
    if (NULL != sf_getstring("flag")) {
	flag = sf_output("flag");
	sf_settype(flag,SF_INT);
	sf_putint(flag,"n3",n[2]);
	sf_putfloat(flag,"d3",d[2]);
	sf_putfloat(flag,"o3",o[2]);

	f = sf_intalloc(nt);
    } else {
	flag = NULL;
	f = NULL;
    }

    /* output DSR characteristic angle */
    if (NULL != sf_getstring("alpha")) {
	alpha = sf_output("alpha");
	sf_putint(alpha,"n3",n[2]);
	sf_putfloat(alpha,"d3",d[2]);
	sf_putfloat(alpha,"o3",o[2]);

	al = sf_floatalloc(nt);
    } else {
	alpha = NULL;
	al = NULL;
    }

    /* allocate memory for output */
    t = sf_floatalloc(nt);

    if (!sf_getbool("causal",&causal)) causal=true;
    /* if y, neglect non-causal branches of DSR */

    /* initialize */
    dsreiko_init(n,o,d,
		 thres,tol,nloop,
		 causal,dp);

    /* compute */
    dsreiko_fastmarch(t,s,f,al);

    /* mirror */
    dsreiko_mirror(t);

    if (flag != NULL) sf_intwrite(f,nt,flag);
    if (alpha != NULL) sf_floatwrite(al,nt,alpha);
    
    /* write output */
    sf_floatwrite(t,nt,out);

    exit(0);
}
