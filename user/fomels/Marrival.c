/* Multiple-arrival interpolation from down-marching. */
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
			     
static sf_eno tfnt, pfnt;
static int it;
static float sx;

static int compfunc(const void *a, const void *b);
static float func_eno(float t);

int main (int argc, char* argv[])
{
    int nt,nx,nz, ig, ix,iz, ng, nw;
    float t, a, b, dz, f, g;
    float *tx, *px, *zx;
    sf_file in, out, place, depth;

    sf_init (argc,argv);
    in = sf_input("in");
    place = sf_input("place");
    out = sf_output("out");

    if (NULL == sf_getstring("depth")) {
	depth = NULL;
    } else {
	depth = sf_input("depth");
    }

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nz)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&dz)) sf_error("No n3= in input");

    if (!sf_getfloat ("sx",&sx)) sx=0.;
    /* source x position */
    if (!sf_getint ("nw",&nw)) nw=3;
    /* interpolation accuracy */

    tx = sf_floatalloc(nt);
    px = sf_floatalloc(nt);
    zx = sf_floatalloc(nt);

    if (NULL == depth) {
	for (it = 0; it < nt; it++) {
	    zx[it] = 0.;
	}
    }

    tfnt = sf_eno_init (nw, nt);
    pfnt = sf_eno_init (nw, nt);

    ng = 0;
    /* iy = nx*nz; */
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    sf_floatread(tx,nt,in);
	    sf_floatread(px,nt,place);
	    if (NULL != depth) sf_floatread(zx,nt,depth);

	    sf_eno_set (tfnt, tx);
	    sf_eno_set (pfnt, px);

	    ig = 0;
	    for (it = 0; it < nt-1; it++) {
		a = px[it]-sx;
		b = px[it+1]-sx;

		if (SF_SIG(a) != SF_SIG(b) && (zx[it] < dz || zx[it+1] < dz)) {
		    t = sf_zero(func_eno,0.,1.,a,b,SF_EPS,false);
		    
		    sf_eno_apply (tfnt,it,t,&f,&g,FUNC);
		    
		    tx[ig] = f;
		    ig++;
		}
	    }        
	    if (ig > ng) ng = ig;
	    if (ig == 0) { /* no arrivals */
		for (it = 0; it < nt; it++) {
		    tx[it] = -1.;
		}
	    } else {
		qsort(tx, ig, sizeof(float), compfunc);
		for (it = ig; it < nt; it++) {
		    tx[it] = -1.;
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
    sf_eno_apply (pfnt,it,t,&f,&g,FUNC);
    return (f-sx);
}

