#include <math.h>

#include <rsf.h>

#include "symes.h"

int main(int argc, char* argv[])
{
    int nx, nz, ix, iz;
    float* trace;
    float dx, dz, x[2];
    sf_file mod;

    sf_init (argc,argv);
    mod = sf_output("out");

    if (!sf_getint("nx",&nx)) nx=400; dx = 1./(nx-1);
    if (!sf_getint("nz",&nz)) nz=800; dz = 2./(nz-1);

    sf_putint   (mod,"n1",nz); 
    sf_putfloat (mod,"d1",dz); 
    sf_putfloat (mod,"o1",0.);
    sf_putint   (mod,"n2",nx); 
    sf_putfloat (mod,"d2",dx); 
    sf_putfloat (mod,"o2",0.);
    sf_setformat (mod,"native_float");

    trace = sf_floatalloc(nz);

    for (ix = 0; ix < nx; ix++) {
	x[1] = ix*dx;
	for (iz = 0; iz < nz; iz++) {
	    x[0] = iz*dz;
	    trace[iz] = 1./sqrtf(symes_vel(NULL,x));
	}
	sf_write(trace,sizeof(float),nz,mod);
    }

    exit (0);
}

