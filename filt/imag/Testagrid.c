#include <math.h>

#include <rsf.h>

#include "agrid.h"

static float x0, dx;

static void mysin (float q, void* v, float* f)
{
  float x;

  x = x0+q*dx;

  f[0] = x;
  f[1] = sinf(x);

  sf_warning("got %g %g\n",f[0],f[1]);
}

int main(int argc, char* argv[])
{
  int maxsplit, nx, ng;
  float **place, **out, min1, max1;
  agrid grd;
  sf_file grid;

  sf_init (argc,argv);
  grid = sf_output("out");
  sf_setformat(grid,"native_float");

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
  
  ng = grid_size(grd);
  out = write_grid(grd);

  sf_putint(grid,"n1",2);
  sf_putint(grid,"n2",ng);

  sf_floatwrite(out[0],ng*2, grid);

  exit(0);
}
