/* Prestack first-arrival traveltime tomography (DSR) */
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

#include "dsrtomo.h"

int main(int argc, char* argv[])
{
    bool adj;
    int dim, i, n[SF_MAX_DIM], it, nt;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], *dt, *dw, *t;
    char key[6];
    sf_file in, out, time;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */

    /* read in time table */
    if (NULL == sf_getstring("time"))
	sf_error("Need time table time=");
    time = sf_input("time");

    /* read dimension */
    dim = sf_filedims(time,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(time,key,d+i)) sf_error("No %s= in input.",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(time,key,o+i)) o[i]=0.;
	nt *= n[i];
    }
    
    /* read input */
    if (adj) {
	dt = sf_floatalloc(nt);
	sf_floatread(dt,nt,in);
    } else {
	dt = sf_floatalloc(nt/n[2]);
	sf_floatread(dt,nt/n[2],in);
    }

    /* output dimension */
    if (adj) {
	sf_putint(out,"n3",1);
    } else {
	sf_putint(out,"n3",n[1]);
	sf_putfloat(out,"d3",d[1]);
	sf_putfloat(out,"o3",o[1]);
    }

    /* output array */
    if (adj) {
	dw = sf_floatalloc(nt/n[2]);
    } else {
	dw = sf_floatalloc(nt);
    }

    /* read stencil time */
    t = sf_floatalloc(nt);
    sf_floatread(t,nt,time);
    sf_fileclose(time);

    /* initialize */
    dsrtomo_init(dim,n,d);

    /* set operator */
    dsrtomo_set(t);

    /* operator */
    if (adj)
	dsrtomo_oper(true,false,n[0]*n[1],n[0]*n[1]*n[1],dw,dt);
    else
	dsrtomo_oper(false,false,n[0]*n[1],n[0]*n[1]*n[1],dw,dt);

    /* write output */
    if (adj) {
	sf_floatwrite(dw,nt/n[2],out);
    } else {
	sf_floatwrite(dw,nt,out);
    }

    exit(0);
}
