/* 1D LUM operator */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include <stdio.h>
#include <math.h>
#include <rsf.h>

#include "lum.h"
#include "median.h"

float lum(float *temp, int nfw, int smnclip, int shnclip)
/*< get a LUM value >*/
{   
  float lsmoother,usmoother,smlow,smupper,shlow,shupper;
  float center,midpoint,output;
  float lum[3];
  center=temp[nfw/2];
  smlow=sf_quantile(smnclip-1,nfw,temp);
  smupper=sf_quantile(nfw-smnclip,nfw,temp);
  shlow=sf_quantile(shnclip-1,nfw,temp);
  shupper=sf_quantile(nfw-shnclip,nfw,temp);
  lum[0]=smlow;
  lum[1]=center;
  lum[2]=shlow;
  lsmoother=medianfilter(lum,3);
  lum[0]=smupper;
  lum[1]=center;
  lum[2]=shupper;
  usmoother=medianfilter(lum,3);
  midpoint=(lsmoother+usmoother)/2.;
  if(center<=midpoint){
    output=lsmoother;
  }else{
    output=usmoother;
  }

  return output;
}

/* 	$Id$	 */


