/* Anisotropic diffusion by regularized inversion. 

Applied in 3D in a slice by slice fashion: set of
2D diffusions for a fixed time sample. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <assert.h>

#include <rsf.h>

#include "anisodiffuse.h"

int main(int argc, char* argv[])
{
    int niter, repeat, n12, n123, n1,n2,n3;
    float *data, *data2, eps, *vx, *vy;
    sf_file in, out, fvx, fvy;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* Get dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2=");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3=");
    n12 = n1*n2;
    n123 = n12*n3;

    /* Allocate and read data */
    data = sf_floatalloc(n123);
    data2 = sf_floatalloc(n123);
    sf_floatread(data,n123,in);

    /* Allocate and read orientation vector */

    fvx = sf_input("vx");
    fvy = sf_input("vy");

    vx = sf_floatalloc(n123);
    sf_floatread(vx,n123,fvx);
    vy = sf_floatalloc(n123);
    sf_floatread(vy,n123,fvy);
    
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of conjugate-gradient iterations */
    if (!sf_getint("repeat",&repeat)) repeat=1;
    /* number of smoothing iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* regularization parameter */

    sf_warning("before init");
    
    anisodiffuse_init(n1,n2,n3 /* data size */, 
		      vx,vy /* parallel to edges vector components */,
		      niter /* number of iterations */,
		      repeat /* number of smoothing iterations */,
		      eps /* regularization */);

    sf_warning("before lop");

    anisodiffuse_lop(n123,n123,data,data2);

    sf_floatwrite(data2,n123,out);

    exit(0);
}
