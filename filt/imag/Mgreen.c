/* Phase-space Green\'s function from down-marching.

Takes < input.rsf place=place.rsf depth=depth.rsf wave=wave.rsf > output.rsf
*/

#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int it, nt,nx,nz, ix,iz, iw, nw;
    float a, dz, z, dw, w0, w1, w2, w, sx, sz, dx;
    float *tx, *px, *zx;
    float complex *trace1, *trace2, c;
    sf_file in, out, place, depth, wave;

    sf_init (argc,argv);
    in = sf_input("in");
    place = sf_input("place");
    depth = sf_input("depth");
    out = sf_output("out");
    wave = sf_input("wave");

    if (SF_COMPLEX != sf_gettype(wave)) sf_error("Need complex wave");
    if (!sf_histint(wave,"n1",&nw)) sf_error("No n1= in wave");
    if (!sf_histfloat(wave,"d1",&dw)) sf_error("No d1= in wave");
    if (!sf_histfloat(wave,"o1",&w0)) sf_error("No o1= in wave");
    if (!sf_getfloat("w1",&w1)) w1=w0;
    /* lowest frequency */
    if (!sf_getfloat("w2",&w2)) w2=w0+(nw-1)*dw;
    /* highest frequency */

    trace1 = sf_complexalloc(nw);
    trace2 = sf_complexalloc(nw);
 
    sf_read(trace1,sizeof(float complex),nw,wave);

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nz)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dz)) sf_error("No d3= in input");
    dz *= (0.5*nz);

    if (!sf_getfloat ("sx",&sx)) sx=0.;
    /* source x position */
    if (!sf_getfloat ("sz",&sz)) sz=0.;
    /* source z position */
    if (!sf_getfloat ("z",&z)) z=0.5;
    /* depth */

    sf_putint(out,"n1",nw);
    sf_putfloat(out,"d1",dw);
    sf_putfloat(out,"o1",w0);
    sf_settype(out,SF_COMPLEX);

    tx = sf_floatalloc(nt);
    px = sf_floatalloc(nt);
    zx = sf_floatalloc(nt);

    for (iz=0; iz<nz; iz++) {
	sf_warning("depth %d of %d",iz+1, nz);
	for (ix=0; ix<nx; ix++) {
	    sf_read(tx,sizeof(float),nt,in);
	    sf_read(px,sizeof(float),nt,place);
	    sf_read(zx,sizeof(float),nt,depth);
	    
	    for (iw = 0; iw < nw; iw++) {
		trace2[iw] = 0.;
		
		w = w0 + iw*dw;
		if (w < w1 || w > w1) continue;
		
		c = trace1[iw];
		w *= 2.*SF_PI;
		
		for (it = 0; it < nt; it++) {
		    if (zx[it] > sz+dz || zx[it] < sz-dz) continue;
		    
		    a = fabsf(px[it]-sx);
		    if (a < dx) {
			trace2[iw] += (1.-a)*c*cexp(I*w*tx[it]);
		    }
		} /* nt */
	    } /* nw */
	    sf_write (trace2,sizeof(float complex),nw,out);
	} /* nx */
    } /* nz */
    
    exit (0);
}

/* 	$Id: Mgreen.c,v 1.3 2003/10/23 03:38:28 fomels Exp $	 */
