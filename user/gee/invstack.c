/* NMO stack by inverse of forward modeling */
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

#include "invstack.h"
#include "imospray.h"

void invstack(int nt, float *model, int nx, const float *gather, 
	      float t0, float x0, 
	      float dt, float dx, float slow, int niter) 
/*< NMO stack by inverse of forward modeling */
{
    imospray_init( slow, x0,dx, t0,dt, nt, nx);
    sf_tinysolver( imospray_lop, sf_cgstep, 
		   nt, nt*nx, model, NULL, gather, niter);
    sf_cgstep_close ();  
    imospray_close ();  /* garbage collection */
}
