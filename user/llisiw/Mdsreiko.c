/* Double square-root eikonal solver */
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
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *s, *t;
    char key[6];
    sf_file in, out;

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
	/* extend the third dimension for output */
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

    /* allocate memory for output */
    t = sf_floatalloc(ns*n[1]);

    /* initialize */
    dsreiko_init(n,o,d);

    /* compute */
    dsreiko_fastmarch(t,s);

    /* write output dimension */
    sf_putint(out,"n3",n[1]);
    sf_putfloat(out,"d3",d[1]);
    sf_putfloat(out,"o3",o[1]);

    /* write output */
    sf_floatwrite(t,ns*n[1],out);

    exit(0);
}
