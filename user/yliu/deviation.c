/* Deviation. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "deviation.h"
#include "mean.h"

float sdeviation(float *temp, int nfw)
/*< calculate a standard deviation >*/
{
    int i;
    float m,sd,*data;
    data = sf_floatalloc(nfw);
    for(i=0;i<nfw;i++){
	data[i]=temp[i];
    }
    m=mean(data,nfw);
    for(i=0;i<nfw;i++){
	data[i]=(data[i]-m)*(data[i]-m);
    }
    sd=mean(data,nfw);
    sd=sqrtf(sd);
    return sd;
}

/* 	$Id$	 */


