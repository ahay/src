/* 2-D synthetic model for multiple-arrival generation.

Takes: > model.rsf

From Dave Hale.
*/

#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nz, nx, ix, iz;
    float *trace, x0, z0, x, z, dx, dz, v0, g, v, a, d;
    sf_file mod;

    sf_init (argc,argv);
    mod = sf_output("out");

    if (!sf_getint("nz",&nz)) nz=601; 
    /* vertical dimension */
    if (!sf_getint("nx",&nx)) nx=401; 
    /* horizontal dimension */

    dz = 4./(nz-1);
    dx = 6./(nx-1);
    
    sf_putint   (mod,"n1",nz); 
    sf_putfloat (mod,"d1",dz); 
    sf_putfloat (mod,"o1",0.);
    sf_putint   (mod,"n2",nx); 
    sf_putfloat (mod,"d2",dx); 
    sf_putfloat (mod,"o2",0.);
    sf_setformat (mod,"native_float");

    if(!sf_getfloat ("x0",&x0)) x0=4.;
    /* anomaly center in x */
    if(!sf_getfloat ("z0",&z0)) z0=1.5;
    /* anomaly center in z */
    if(!sf_getfloat ("v0",&v0)) v0=1.5;
    /* surface velocity */
    if(!sf_getfloat ("g",&g)) g=0.6;
    /* velocity gradient */
    if(!sf_getfloat ("a",&a)) a=1.;
    /* anomaly magnitude */
    if(!sf_getfloat ("d",&d)) d=1.;
    /* anomaly radius */
    d *= d;

    trace = sf_floatalloc(nz);

    for (ix=0; ix <  nx; ix++) {
	x = ix*dx - x0;
	x *= x;
	for (iz=0; iz < nz; iz++) {
	    z = iz*dz;
	    v = v0 + g*z;
	    z = z - z0;
	    z = z*z + x;
	    trace[iz] = v + a*exp(-z/d);
	}
	sf_floatwrite (trace,nz,mod);
    }

    exit (0);
}

/* 	$Id: Mhale.c,v 1.2 2004/04/19 21:51:16 fomels Exp $	 */


