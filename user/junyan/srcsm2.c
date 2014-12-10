/* All-pass plane-wave destruction filter coefficients */
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

#include "srcsm2.h"


static float **wt;

void srcsm_init(float dz, float dx/*grid size*/)
/*< initialization >*/
{

    float dx2, dz2, R0, rxxz, rxzz;
    dx2 = dx*dx;
    dz2 = dz*dz;
    R0 = sqrtf(dx2+dz2);
    rxxz = sqrtf(4.0*dx2+dz2);
    rxzz = sqrtf(4.0*dz2+dx2);
    wt =  sf_floatalloc2(3,2);
    wt[0][0] = -0.5*dz2/(R0*R0)+1.0;
    if(2.0*dz < R0) { 
      wt[0][1]= -2.0*dz2/(R0*R0)+1.0;
    } else {
      wt[0][1]= 2.0*(dz-R0)*(dz-R0)/(R0*R0);
    }
    wt[0][2]= 0.5*(rxzz-2.0*R0)*(rxzz-2.0*R0)/(R0*R0);
    wt[1][0] = -0.5*dx2/(R0*R0)+1.0;
    if(2.0*dx < R0) { 
      wt[1][1] = -2.0*dx2/(R0*R0)+1.0;
    } else {
      wt[1][1] = 2.0*(dx-R0)*(dx-R0)/(R0*R0);
    }
    wt[1][2]= 0.5*(rxxz-2.0*R0)*(rxxz-2.0*R0)/(R0*R0);



}
   

void srcsm_close(void)
/*< free memory allocation>*/
{
    free(*wt);
    free(wt);

}



void source_smooth(float **source /*source matrix*/, 
                   int iz /*source index*/, 
                   int ix /*source index*/,
                   float val/*source value*/)
/*< smooth source input >*/
{
   
    source[iz-1][ix]   += val*wt[0][0]; 
    source[iz+1][ix]   += val*wt[0][0]; 

    source[iz+2][ix]   += val*wt[0][1];
    source[iz-2][ix]   += val*wt[0][1];

    source[iz+2][ix+1]   += val*wt[0][2]; 
    source[iz+2][ix-1]   += val*wt[0][2]; 
    source[iz-2][ix-1]   += val*wt[0][2]; 
    source[iz-2][ix+1]   += val*wt[0][2]; 

    source[iz-1][ix-1] += val*0.5; 
    source[iz-1][ix+1] += val*0.5; 
    source[iz+1][ix-1] += val*0.5; 
    source[iz+1][ix+1] += val*0.5; 

    source[iz][ix-1]   += val*wt[1][0];
    source[iz][ix+1]   += val*wt[1][0];

    source[iz][ix+2]   += val*wt[1][1];
    source[iz][ix-2]   += val*wt[1][1];

    
    source[iz+1][ix+2]   += val*wt[1][2]; 
    source[iz+1][ix-2]   += val*wt[1][2]; 
    source[iz-1][ix-2]   += val*wt[1][2]; 
    source[iz-1][ix+2]   += val*wt[1][2]; 

}

void source_smooth1(float dz /*grid size*/,
                    float dx /*grid size*/, 
                    float alpha /*decay parameter*/, 
                    float **source /*source matrix*/, 
                    int iz /*source index*/, 
                    int ix /*source index*/,
                    float val/*source value*/)
/*< smooth source input >*/
{
    float r, dx2, dz2, l;
    r = SF_MIN(dx,dz);
    dx2 = dx*dx;
    dz2 = dz*dz;
    l = sqrtf(dx2+dz2);
    source[iz-1][ix]   += val*exp(alpha*(dz/r)); 
    source[iz+1][ix]   += val*exp(alpha*(dz/r)); 
    source[iz][ix-1]   += val*exp(alpha*(dx/r)); 
    source[iz][ix+1]   += val*exp(alpha*(dx/r)); 
    source[iz-1][ix-1] += val*exp(alpha*(l/r)); 
    source[iz-1][ix+1] += val*exp(alpha*(l/r)); 
    source[iz+1][ix-1] += val*exp(alpha*(l/r)); 
    source[iz+1][ix+1] += val*exp(alpha*(l/r)); 

}
