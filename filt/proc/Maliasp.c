/* Aliasing test. 

Takes: > aliased.rsf
*/

#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt,nx,it,ix,ix0;
    float *wave, *data, cycles;
    sf_file out;
    
    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");

    if (!sf_getint("n1",&nt)) nt=600;
    if (!sf_getint("n2",&nx)) nx=24;
    /* dimensions */

    if (!sf_getfloat("cycles",&cycles)) cycles=10.;
    /* wave frequency */

    if (!sf_getint("ix0",&ix0)) ix0=0; 
    /* central trace */
    /* try ix0=2 */

    sf_putint(out,"n1",nt);
    sf_putint(out,"n2",nx);
    sf_putint(out,"d1",1.);
    sf_putint(out,"d2",1.);

    data = sf_floatalloc(nt);
    wave = sf_floatalloc(nt);

    for (it=0; it < nt; it++) {
	data[it]=0.;
	wave[it] = sinf(2.*SF_PI*it*cycles/nt) * 
	    expf(- 3.*(it+1.)/nt);
    }

    for (ix=0; ix < nx; ix++) {
	it = nt*sqrtf(1. + 0.01*(ix-ix0)*(ix-ix0))/3.;
	if (it > nt) it=nt;
	if (it > 0)  sf_write(data,sizeof(float),it,out);
	if (nt > it) sf_write(wave,sizeof(float),nt-it,out);
    }
 
    exit(0);
}
