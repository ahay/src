/* Multiple arrivals by depth marching.

Takes: < velocity.rsf > arrivals.rsf place=place.rsf depth=depth.rsf

*/

#include <math.h>

#include <rsf.h>

#include "ztrace.h"

int main(int argc, char* argv[])
{
    bool vel;
    int nz,nx, iz, na, nt, nax, ix, order, is, iorder, k;
    float **slow, *slice[NS];
    float dz,dx,da,x0,z0,a0;
    sf_file in, time, place, depth;

    sf_init(argc,argv);
    in = sf_input("in");
    time = sf_output("out");
    place = sf_output("place");
    depth = sf_output("depth");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
  
    if (!sf_histfloat(in,"o1",&z0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_getint("nt",&nt)) nt=nx*nz;
    /* ray length bound */

    sf_putint(time,"n3",nz);
    sf_putfloat(time,"d3",dz);
    sf_putfloat(time,"o3",z0);

    sf_putint(time,"n2",nx);
    sf_putfloat(time,"d2",dx);
    sf_putfloat(time,"o2",x0);

    if (!sf_getint("na",&na)) na=362;
    /* number of angles */
    if (!sf_getfloat("da",&da)) da=0.5;
    /* angle increment (in degrees) */
    if (!sf_getfloat("a0",&a0)) a0=-90.;
    /* starting angle (in degrees) */

    sf_putint(time,"n1",na);
    sf_putfloat(time,"d1",da);
    sf_putfloat(time,"o1",a0);

    sf_putint(place,"n1",na);
    sf_putfloat(place,"d1",da);
    sf_putfloat(place,"o1",a0);

    sf_putint(depth,"n1",na);
    sf_putfloat(depth,"d1",da);
    sf_putfloat(depth,"o1",a0);

    da *= (SF_PI/180.);
    a0 *= (SF_PI/180.);

    nax = na*nx;

    sf_putint(place,"n3",nz);
    sf_putfloat(place,"d3",dz);
    sf_putfloat(place,"o3",z0);

    sf_putint(place,"n2",nx);
    sf_putfloat(place,"d2",dx);
    sf_putfloat(place,"o2",x0);

    sf_putint(depth,"n3",nz);
    sf_putfloat(depth,"d3",dz);
    sf_putfloat(depth,"o3",z0);

    sf_putint(depth,"n2",nx);
    sf_putfloat(depth,"d2",dx);
    sf_putfloat(depth,"o2",x0);
    
    /* additional parameters */
    if(!sf_getbool("vel",&vel)) vel=true;
    /* y, input is velocity; n, slowness */
    if(!sf_getint("order",&order)) order=3;
    /* interpolation accuracy for velocity */
    if(!sf_getint("iorder",&iorder)) iorder=4;
    /* interpolation accuracy for grid */

    slow  = sf_floatalloc2(nz,nx);
    sf_read(slow[0],sizeof(float),nz*nx,in);
    if (vel) { /* convert to slowness */
	for(ix = 0; ix < nx; ix++){
	    for (iz = 0; iz < nz; iz++) {
		slow[ix][iz] = 1./slow[ix][iz];
	    }
	}
    }
    
    for (is = 0; is < NS; is++) {
	slice[is] = sf_floatalloc(nax);
	for (k = 0; k < nax; k++) {
	    slice[is][k] = 0.;
	}
    }

    ztrace_init (order, iorder, nx, nz, na, nt,
                 dx, dz, da, x0, z0, a0, slow, slice);

    for (iz = 0; iz < nz; iz++) {
	sf_warning("depth %d of %d", iz+1, nz);

	ztrace_step (iz);

	sf_write (slice[0],sizeof(float),nax,time);
	sf_write (slice[1],sizeof(float),nax,place);
	sf_write (slice[2],sizeof(float),nax,depth);
    }

    exit (0);
}

/* 	$Id: Mztrace.c,v 1.1 2003/10/21 15:09:08 fomels Exp $	 */
