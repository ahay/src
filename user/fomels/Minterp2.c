/* Multiple-arrival interpolation. */
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
static float sx, sz;

static int compfunc(const void *a, const void *b);
static float func_eno(float t);

int main (int argc, char* argv[])
{
    int four,nt,nx,nz, ig, ix,iz, ng, nw, nt2,st, plane, **siz;
    float t, a, b, f, g, dz, xztp[5];
    float *tx, *px, *zx;
    sf_file in, out, size, grid;

    sf_init (argc,argv);
    in = sf_input("in");
  
    if (!sf_histint(in,"n1",&four)) sf_error("No n1= in input");
    if (!sf_histint(in,"n3",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n4",&nz)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d4",&dz)) sf_error("No d3= in input");

    if (!sf_getfloat ("sx",&sx)) sx=0.;
    if (!sf_getfloat ("sz",&sz)) sz=0.;
    /* Shot coordinates */
    if (!sf_getint ("nw",&nw)) nw=4;
    /* Interpolation accuracy */
    if (!sf_getint ("plane",&plane)) plane=0;
    /* 0: point-source, 4: plane-wave */

    size = sf_input("size");
    siz = sf_intalloc2 (nx,nz);
    sf_intread(siz[0],nx*nz,size);
    sf_fileclose(size);

    nt = 0;
    for (iz=0; iz<nz; iz++) {
	for (ix=0; ix<nx; ix++) {
	    st = siz[iz][ix];
	    if (nt < st) nt=st;
	}
    }
    sf_warning("maxsize is %d",nt);

    out = sf_output("out");
    sf_putint(out,"n1",nt);
    sf_putint(out,"n2",1);

    grid = sf_input("grid");

    tx = sf_floatalloc(nt);
    px = sf_floatalloc(nt);
    zx = sf_floatalloc(nt);

    ng = 0;
    for (iz=0; iz<nz; iz++) {
	sf_warning("depth %d of %d;",iz+1, nz);
    
	for (ix=0; ix<nx; ix++) {
	    nt2 = siz[iz][ix];

	    for (it=0; it < nt2; it++) {
		sf_floatread(xztp, four, grid);
		tx[it] = xztp[2];
		px[it] = xztp[plane];
		zx[it] = xztp[1];
	    }

	    tfnt = sf_eno_init (nw, nt2);
	    pfnt = sf_eno_init (nw, nt2);

	    sf_eno_set (tfnt, tx);
	    sf_eno_set (pfnt, px);
 
	    ig = 0;
	    for (it = 0; it < nt2-1; it++) {
		if (zx[it] > sz+dz || zx[it+1] > sz+dz) continue;

		a = px[it]-sx;
		b = px[it+1]-sx;
	
		if ((a <= 0. && b > 0.) ||
		    (a >= 0. && b < 0.)) {
	  
		    t = sf_zero(func_eno,0.,1.,a,b,1.e-3,false);
		    sf_eno_apply (tfnt,it,t,&f,&g,FUNC);
	  
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
		for (it = ig; it < nt; it++) {
		    tx[it] = -1.;
		}
	    }
	    sf_floatwrite (tx,nt,out);
      
	    sf_eno_close (tfnt);
	    sf_eno_close (pfnt);
	}
    }
    sf_warning(".");
  
    sf_warning("number of branches = %d", ng);

    sf_close();
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

/* 	$Id: Minterp2.c 4084 2009-01-22 21:39:08Z sfomel $	 */
