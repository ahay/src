/* Phase-space Green\'s function from down-marching.

Takes < input.rsf place=place.rsf depth=depth.rsf wave=wave.rsf > output.rsf
*/

#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int it, nt,nx,nz, ix,iz, iw, nw;
    float a, dz, z, dw, w0, w1, w2, w, sx, sz, dx;
    float *tx, *px, *zx, *trace2;
    float complex *trace1, c;
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
    trace2 = sf_floatalloc(nx);
 
    sf_read(trace1,sizeof(float complex),nw,wave);

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nz)) sf_error("No n3= in input");

    if (!sf_getfloat("dx",&dx)) sf_error("Need dx= ");
    /* location error */

    if (!sf_histfloat(in,"d3",&dz)) sf_error("No d3= in input");
    dz *= (0.5*nz);

    if (!sf_getfloat ("sx",&sx)) sx=0.;
    /* source x position */
    if (!sf_getfloat ("sz",&sz)) sz=0.;
    /* source z position */
    if (!sf_getfloat ("z",&z)) z=0.5;
    /* depth */

    sf_putint(out,"n1",1);
/*    
      sf_putfloat(out,"d1",dw);
      sf_putfloat(out,"o1",w0); 
      sf_settype(out,SF_COMPLEX); 
*/

    tx = sf_floatalloc(nt);
    px = sf_floatalloc(nt);
    zx = sf_floatalloc(nt);

    for (iz=0; iz<nz; iz++) {
	sf_warning("depth %d of %d",iz+1, nz);
	for (ix=0; ix<nx; ix++) {
	    sf_read(tx,sizeof(float),nt,in);
	    sf_read(px,sizeof(float),nt,place);
	    sf_read(zx,sizeof(float),nt,depth);
	    
	    trace2[ix] = 0.;

	    for (iw = 0; iw < nw; iw++) {
		w = w0 + iw*dw;
		if (w < w1 || w > w2) continue;
		
		c = trace1[iw];
		if (fabsf(w) < dw) c *= 2.;
		w *= 2.*SF_PI;
		
		for (it = 0; it < nt; it++) {
		    if (zx[it] > sz+dz || zx[it] < sz-dz) continue;
		    
		    a = fabsf(px[it]-sx);
		    if (a < dx) {
			trace2[ix] += crealf((1.-a/dx)*c*cexpf(-I*w*tx[it]));
		    }
		} /* nt */
	    } /* nw */
	} /* nx */
	sf_write (trace2,sizeof(float),nx,out);
    } /* nz */
    
    exit (0);
}

/* 	$Id: Mgreen.c,v 1.4 2003/10/24 14:57:58 fomels Exp $	 */
