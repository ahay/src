/* Given a PEF, finds bad data and restores it. */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

#include "mis2.h"

void fixbad (int niter    /* number of iterations */, 
	     sf_filter aa /* PEF */, 
	     int ny       /* data size */,
	     float *yy    /* in - data, out - deburst */) 
/*< find bad data and restore it >*/
{
    int iy;
    float *rr, *rabs, rbar;
    bool *known;

    rr = sf_floatalloc(ny);
    rabs = sf_floatalloc(ny);
    known = sf_boolalloc(ny);

    sf_helicon_init(aa);
    sf_helicon_lop (false,false,ny,ny,yy,rr); 
    for (iy=0; iy < ny; iy++) 
	rabs[iy] = fabsf (rr[iy]);
    rbar = sf_quantile(ny/2,ny,rabs);
    for (iy=0; iy < ny; iy++) 
	known[iy] = (bool) ((yy[iy] > 0.) && (fabsf(rr[iy]) < 4. * rbar));
    mis2 (niter, ny, yy, aa, known, 0., true);

    free(rr);
    free(rabs);
    free(known);
}

