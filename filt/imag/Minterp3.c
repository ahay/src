/* Multiple-arrival interpolation from down-marching.
*/
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

#include "eno.h"
#include "fzero.h"

static eno tfnt, pfnt, zfnt;
static int it;
static float sx, sz;

static int compfunc(const void *a, const void *b);
static float func_eno(float t);

int main (int argc, char* argv[])
{
    int nt,nx,nz, ig, ix,iy,iz, ng, nw;
    float t, a, b, dt, f, g, dz, z;
    float *tx, *px, *zx;
    sf_file in, out, place, depth;

    sf_init (argc,argv);
    in = sf_input("in");
    place = sf_input("place");
    depth = sf_input("depth");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nz)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dz)) sf_error("No d3= in input");
    dz *= (0.5*nz);

    if (!sf_getfloat ("sx",&sx)) sx=0.;
    /* source x position */
    if (!sf_getfloat ("sz",&sz)) sz=0.;
    /* source z position */
    if (!sf_getfloat ("dt",&dt)) dt=2.e-3;
    /* time uncertainty */
    if (!sf_getint ("nw",&nw)) nw=4;
    /* interpolation accuracy */
    if (!sf_getfloat ("z",&z)) z=0.5;
    /* depth */

    tx = sf_floatalloc(nt);
    px = sf_floatalloc(nt);
    zx = sf_floatalloc(nt);

    tfnt = eno_init (nw, nt);
    pfnt = eno_init (nw, nt);
    zfnt = eno_init (nw, nt);

    ng = 0;
    iy = nx*nz;
    for (iz=0; iz<nz; iz++) {
	sf_warning("depth %d of %d",iz+1, nz);

	for (ix=0; ix<nx; ix++) {
	    sf_floatread(tx,nt,in);
	    sf_floatread(px,nt,place);
	    sf_floatread(zx,nt,depth);

	    eno_set (tfnt, tx);
	    eno_set (pfnt, px);
	    eno_set (zfnt, zx);

	    ig = 0;
	    for (it = 0; it < nt-1; it++) {
		if ((zx[it] > sz+dz && zx[it+1] > sz+dz) ||
		    (zx[it] < sz-dz && zx[it+1] < sz-dz)) continue;

		a = px[it]-sx;
		b = px[it+1]-sx;

		if ((a <= 0. && b > 0.) ||
		    (a >= 0. && b < 0.)) {
	    
/*	    fprintf(stderr,"%d %d %d %f %f\n", it, ix, iz, a, b); */

		    t = fzero(func_eno,0.,1.,a,b,1.e-3,false);

		    eno_apply (zfnt,it,t,&f,&g,FUNC);
		    if (f > sz + dz || 
			f < sz - dz) continue;

		    eno_apply (tfnt,it,t,&f,&g,FUNC);
	  
		    if (ix + nx*iz == iy) {
			if (fabs (tx[ig-1]-f) < dt) continue;
		    } else {
			iy = ix + nx*iz;
		    }
	  
		    tx[ig] = f;
		    ig++;
		}
	    }        
	    if (ig > ng) ng = ig;
	    if (ig == 0) {
		for (it = 0; it < nt; it++) {
		    tx[it] = -1.;
		}
	    } else {
		qsort(tx, ig, sizeof(float), compfunc);
		if (ig > 1) {
		    if (ig < nt) tx[ig] = tx[ig-1];
		    for (it = ig+1; it < nt; it++) {
			tx[it] = -1.;
		    }
		} else {
		    for (it = ig; it < nt; it++) {
			tx[it] = -1.;
		    }
		}
	    }
	    sf_floatwrite (tx,nt,out);
	}
    }
  
    sf_warning("number of branches = %d", ng);

    exit (0);
}

static int compfunc(const void *a, const void *b)
{
    float aa, bb;

    aa = *(float *)a;
    bb = *(float *)b;

    if (aa <  bb) return (-1);
    if (aa == bb) return 0;
    return 1;
}


static float func_eno(float t)
{
    float f, g;
    eno_apply (pfnt,it,t,&f,&g,FUNC);
    return (f-sx);
}

/* 	$Id: Minterp3.c,v 1.5 2004/07/02 11:54:20 fomels Exp $	 */
