/* Generic least-squares optimization */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <rsf.h>

void inverts(sf_operator oper /* linear operator */, 
	    int niter        /* number of iterations */,
	    int miter        /* memory */, 
            int liter        /* number of linear iterations */,
	    int nx, int ny   /* model and data size */,
	    float *x         /* model */, 
	    const float *y   /* data */, 
	    float *error     /* least-squares error */,
            bool verb        /* verbosity flag */,
            float eps        /* epsilon */)
/*< Iterative least-squares optimization >*/
{
    int iter, ix;
    float *w, *p;

    sf_cdstep_init();

    w = sf_floatalloc(nx);
    p = sf_floatalloc(nx);

    for (ix=0; ix < nx; ++ix) w[ix] = 1.f;

    for (iter = 0; iter < niter; ++iter) {
        sf_solver_prec (oper, sf_cgstep, sf_copy_lop, nx, nx, ny, x, y, 
                        liter, eps, "err",error, "verb", verb, "mwt", w, "xp", p, "end");

        for (ix=0; ix < nx; ++ix) w[ix] = fabsf(p[ix]);
    }
    sf_cdstep_close();

}
