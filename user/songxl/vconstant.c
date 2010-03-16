/* calculate v0 */
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

#include "vconstant.h"
   
float v0_cal(float *v /* v[nxb] */ ,
             int nxb /* size of v */ , 
             int v0kind /* the kind of v0, 0:rms 1:mean 2:minimum 3:maximum*/)
/*< calculate v0 >*/
{
  int ix;
  float v0=0.0;
  switch(v0kind){
        case(0):/* v0 RMS velocity*/
            for (ix=0; ix < nxb; ix++) v0 += v[ix]*v[ix];
            v0 /= (float)nxb;
            v0 = sqrt(v0);
        break;
        case(1):/* v0 mean velocity*/
            for (ix=0; ix < nxb; ix++) v0 += v[ix];
            v0 /= (float)nxb;
        break;
        case(2):/* v0 minimum velocity*/
            v0 = v[0];
            for (ix=1; ix < nxb; ix++) 
                if(v[ix]<v0) v0 = v[ix];
        break;
        case(3):/* v0 maximum velocity*/
            v0 = v[0];
            for (ix=1; ix < nxb; ix++) 
                if(v[ix]>v0) v0 = v[ix];
        break;
  }
  return(v0);
}
