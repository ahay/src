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

void invert(sf_operator oper /* linear operator */, 
	    int niter        /* number of iterations */,
	    int miter        /* memory */, 
	    int nx, int ny   /* model and data size */,
	    float *x         /* model */,
            float *x0        /* inital model */, 
	    const float *y   /* data */, 
	    float *error     /* least-squares error */)
/*< Iterative least-squares optimization >*/
{
    int iy, iter;
    float norm;

    sf_cdstep_init();
    sf_solver(oper,sf_cdstep,nx,ny,x,y,niter,
	      "x0",x0,"nmem",0,"nfreq",miter,"err",error,"verb",true,"end");
    sf_cdstep_close();

    norm = 0.;
    for (iy=0; iy < ny; iy++) {
	norm += y[iy]*y[iy];
    }
    for (iter=0; iter < niter; iter++) {
	error[iter] /= norm;
    }
}
