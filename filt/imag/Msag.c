#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt, nx, it,ix;
    float t,t0,dt, x,x0,dx, z; 
    float tmax,xmax,delx, dxdz,v0,alpha;
    float** datr;
    sf_file sag;

    sf_init (argc, argv);
    sag = sf_output("out");

    if (!sf_getint ("nt",&nt)) nt = 200;
    if (!sf_getint ("nx",&nx)) nx = 200;

    if (!sf_getfloat ("tmax",&tmax)) tmax = 4.;
    if (!sf_getfloat ("xmax",&xmax)) xmax = 4.;
    if (!sf_getfloat ("delx",&delx)) delx = .5;
    if (!sf_getfloat ("dxdz",&dxdz)) dxdz = 1.;
    if (!sf_getfloat ("v0",&v0)) v0 = 1.5;
    if (!sf_getfloat ("alpha",&alpha)) alpha = 0.;

    t0 = 0;		
    x0 = 0.;
    dt = tmax / nt;
    dx = xmax / nx;

    sf_setformat(sag,"native_float");    
    sf_putint (sag,"n1",nt); 
    sf_putfloat (sag,"o1",t0); 
    sf_putfloat (sag,"d1",dt);
    sf_putint (sag,"n2",nx); 
    sf_putfloat (sag,"o2",x0); 
    sf_putfloat (sag,"d2",dx);
    sf_putstring (sag,"label1","Pseudo-depth (s)");
    sf_putstring (sag,"label2","Lateral (km)");

    datr = sf_floatalloc2(nt,nx);
    for (ix=0; ix < nx; ix++) {
	for (it=0; it < nt; it++) {
	    datr[ix][it] = 0.;
	}
    }

    for (x = delx/2; x <= xmax; x+= delx) {
	z = x / dxdz;
	if( alpha != 0.) {
	    t = 2. * log( 1. + alpha * z / v0) / alpha;
	} else {		 
	    t = 2. * z / v0;
	}
	it = 0.5 + t/dt;
	ix = 0.5 + x/dx;
	if(ix < nx && it < nt) datr[ix][it] += 1.;
    }

    sf_write (datr[0],sizeof(float),nt*nx,sag);

    exit (0);
}

/*	Notes about theory
 *	v = v0 + alpha * z
 *	t = 2 * \int_0^z dz/(v0+alpha*z)
 *	t = 2 * log( 1 + alpha*z/v0 ) / alpha
 *	exp(alpha*tmax/2.) = 1 + alpha * zmax/ v0
 *	v0 * (exp( alpha * tmax/2.) - 1) /alpha =  zmax = dz * (nz+1)
 */
