/* Multiple-arrival interpolation (yet another).

Takes: < input.rsf > output.rsf

*/

#include <math.h>

#include <rsf.h>

#include "eno.h"
#include "fzero.h"

static eno tfnt, pfnt;
static int it;
static float sx, sz;

static int compfunc(const void *a, const void *b);
static float func_eno(float t);

int main (int argc, char* argv[])
{
    int nt,nx,nz, ig, ix,iz, ng, nw, four;
    float t, a, b, f, g, dz;
    float *tx, *px, *zx, *xztp;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&four)) sf_error ("No n1= in input");
    if (!sf_histint(in,"n2",&nt)) sf_error ("No n2= in input");
    if (!sf_histint(in,"n3",&nx)) sf_error ("No n3= in input");
    if (!sf_histint(in,"n4",&nz)) sf_error ("No n4= in input");
    if (!sf_histfloat(in,"d4",&dz)) sf_error ("No d4= in input");
    
    if (!sf_getfloat ("sx",&sx)) sx=0.;
    if (!sf_getfloat ("sz",&sz)) sz=0.;
    /* Shot coordinate */
    if (!sf_getint ("nw",&nw)) nw=4;
    /* Interpolation accuracy */

    sf_putint (out,"n1",1);

    tx = sf_floatalloc(nt);
    px = sf_floatalloc(nt);
    zx = sf_floatalloc(nt);
    xztp = sf_floatalloc(four);
    
    tfnt = eno_init (nw, nt);
    pfnt = eno_init (nw, nt);

    ng = 0;
    for (iz=0; iz<nz; iz++) {
	sf_warning("depth %d of %d",iz+1, nz);
	for (ix=0; ix<nx; ix++) {	    
	    for (it=0; it < nt; it++) {
		sf_read(xztp,sizeof(float),four,in);
		tx[it] = xztp[2];
		px[it] = xztp[0];
		zx[it] = xztp[1];
	    }

	    eno_set (tfnt, tx);
	    eno_set (pfnt, px);
 
	    ig = 0;
	    for (it = 0; it < nt-1; it++) {
		if (zx[it] > sz+dz || zx[it+1] > sz+dz) continue;

		a = px[it]-sx;
		b = px[it+1]-sx;

		if ((a <= 0. && b > 0.) ||
		    (a >= 0. && b < 0.)) {

		    t = fzero(func_eno,0.,1.,a,b,1.e-3,false);
		    eno_apply (tfnt,it,t,&f,&g,FUNC);
	  
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
	    sf_write (tx,sizeof(float),nt,out);
	}
    }
    
    eno_close (tfnt);
    eno_close (pfnt);
  
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
    eno_apply (pfnt,it,t,&f,&g,FUNC);
    return (f-sx);
}

/* 	$Id: Minterpt.c,v 1.5 2004/03/22 05:43:24 fomels Exp $	 */
