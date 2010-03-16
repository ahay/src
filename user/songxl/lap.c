/* 1-D second-order derivative*/
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
#include <rsf.h>
#include <float.h>
#include <stdio.h>
#include "lap.h"

static int nxb;

void lap1_init(int n /* */)
/*< initialize >*/
{
   nxb = n;
}



   
void lap1(float *uin /* input */, 
          float *uout /*output*/, 
          int order /*Laplacian order*/)
/*< 1-D second-order derivative >*/
{
  int ix;
  switch(order){
        case(2):
	uout[0] = uin[1]-2*uin[0];
	for (ix=1; ix < nxb-1; ix++) {
	    uout[ix] = uin[ix+1]+uin[ix-1]-2.0*uin[ix]; 
	}
	uout[nxb-1] = uin[nxb-2]-2.0*uin[nxb-1]; 
        break;
        case(4):
	uout[0] = (-uin[2]+16.0*uin[1]-30.0*uin[0])/12.0; 
	uout[1] = (-uin[3]+16.0*uin[2]+16.0*uin[0]-30.0*uin[1])/12.0; 
	for (ix=2; ix < nxb-2; ix++) {
	    uout[ix] = (-uin[ix-2]-uin[ix+2]+16.0*uin[ix-1]+16.0*uin[ix+1]-30.0*uin[ix])/12.0; 
	}
	uout[nxb-1] = (-uin[nxb-3]+16.0*uin[nxb-2]-30.0*uin[nxb-1])/12.0; 
	uout[nxb-2] = (-uin[nxb-4]+16.0*uin[nxb-3]+16.0*uin[nxb-1]-30.0*uin[nxb-2])/12.0; 
        break;
  }
}
