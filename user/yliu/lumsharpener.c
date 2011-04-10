/* 1D LUM sharpener operator */
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

#include "lumsharpener.h"

float lumsharpener(float *temp, int nfw, int nclip)
/*< get a LUM sharpener value >*/
{   
    float sharpener,low,upper,center,midpoint;
    center=temp[nfw/2];
    low=sf_quantile(nclip-1,nfw,temp);
    upper=sf_quantile(nfw-nclip,nfw,temp);
    midpoint=(low+upper)/2.;
    if(center>low && center<upper){
      if(center<=midpoint){
	sharpener=low;
      } else {
	sharpener=upper;
      }
    } else {
      sharpener=center;
    }
    return sharpener;
}
 
/* 	$Id$	 */


