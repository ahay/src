/* A "cup" synthetic model.

Takes: > cup.rsf
*/
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nz, nx, ix, iz;
    float d, x0, *trace;
    sf_file cup;
  
    sf_init (argc,argv);
    cup = sf_output("out");
    sf_setformat(cup,"native_float");
    
    if (!sf_getint("n1",&nz)) sf_error("Need n1=");
    /* vertical axis */
    sf_putint(cup,"n1",nz);
    if (!sf_getint("n2",&nx)) sf_error("Need n2=");
    /* horizontal axis */
    sf_putint(cup,"n2",nx);
    sf_putint(cup,"n3",1);
    
    if (!sf_getfloat("d1",&d)) sf_error("Need d1=");
    /* vertical sampling */
    sf_putfloat(cup,"d1",d);
    if (!sf_getfloat("d2",&d)) sf_error("Need d2=");
    /* horizontal sampling */
    sf_putfloat(cup,"d2",d);
    if (!sf_getfloat("o1",&d)) sf_error("Need o1=");
    /* vertical origin */
    sf_putfloat(cup,"o1",d);
    if (!sf_getfloat("o2",&d)) sf_error("Need o2=");
    /* horizontal origin */
    sf_putfloat(cup,"o2",d);

    trace = sf_floatalloc(nz);

    x0 = (nx/2 + 1) * 0.5;

    for (iz=0; iz < nz; iz++) {
	trace[iz] = 0.;
    }
    iz=0;

    for (ix=0; ix < nx; ix++) {
	trace[iz] = 0.;
	iz = nz/9.0*(3.+cosf((ix/x0-1.0)*SF_PI));
	trace[iz] = 1.;
	sf_write (trace,sizeof(float),nz,cup);
    }

    exit(0);
}

