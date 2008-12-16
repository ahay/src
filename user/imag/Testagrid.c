#include <math.h>

#include <rsf.h>

#include "agrid.h"

static float x0, dx;

static void mysin (float q, void* v, float* f)
{
  float x;

  x = x0+q*dx;

  f[0] = x;
  /* f[1] = sinf(x); */
  f[1] = x*x;
}

int main(int argc, char* argv[])
{
    int maxsplit, i, ix, nx, ng, ig;
    float **place, **out, min1, max1, x, f[2];
    sf_complex *c;
    agrid grd;
    sf_file grid;
    
    sf_init (argc,argv);
    grid = sf_output("out");
    sf_setformat(grid,"native_complex");
    
    if (!sf_getint("nx",&nx)) nx=5;
    /* number of points */
    if (!sf_getfloat("dx",&dx)) dx=SF_PI/(nx-1);
    /* increment */
    if (!sf_getfloat("x0",&x0)) x0=0.;
    /* origin */
    
    if (!sf_getint("maxsplit",&maxsplit)) maxsplit=10;
    /* maximum splitting for adaptive grid */
    
    if (!sf_getfloat("min",&min1)) min1=0.1;
    /* parameters for adaptive grid */
    if (!sf_getfloat("max",&max1)) max1=0.2;
    
    place = sf_floatalloc2(2,nx);
    
    grd = agrid_init (nx, 2, maxsplit);
    agrid_set (grd,place);
    
    fill_grid(grd,0.,2.*dx,min1,max1,NULL,mysin);
    
    for (i=0; i < 4; i++) {
	x = i/dx;
	ix = floorf(x);
	x -= ix;
	agrid_interp(grd,ix,x,f); 
	sf_warning("%g^2 = %g^2 = %g",ix*dx+x,f[0],f[1]);
    }
    
    ng = grid_size(grd);
    out = write_grid(grd);
    
    sf_putint(grid,"n1",ng);
    
    c = sf_complexalloc(ng);
    
    for (ig=0; ig < ng; ig++) {
	c[ig] = sf_cmplx(out[ig][0],out[ig][1]);
    }
    
    sf_complexwrite(c, ng, grid);
    
    exit(0);
}
