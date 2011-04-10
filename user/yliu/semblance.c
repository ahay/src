/* Calculate semblance. */
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

#include "semblance.h"
#include "mean.h"

float semblance(float *temp, int nfw1, int nfw2)
/*< calculate a semblance in a rectangle data >*/
{
    int i,j;
    float sem,*data,*tdata,numerator,denominator;
    data = sf_floatalloc(nfw1*nfw2);
    tdata = sf_floatalloc(nfw1);
    for(i=0;i<nfw1*nfw2;i++){
	data[i]=temp[i];
    }
    for(i=0;i<nfw1;i++){
	tdata[i] = 0.;
	for(j=0;j<nfw2;j++){
	    tdata[i] += data[j*nfw1+i];
	}
    }
    numerator = 0.;
    for(i=0;i<nfw1;i++){
	numerator += tdata[i]*tdata[i];
    }
    for(i=0;i<nfw1;i++){
	tdata[i] = 0.;
	for(j=0;j<nfw2;j++){
	    tdata[i] += data[j*nfw1+i]*data[j*nfw1+i];
	}
    }
    denominator = 0.;
    for(i=0;i<nfw1;i++){
	denominator += tdata[i];
    }
    denominator *= nfw2;

    sem = numerator/(denominator+FLT_EPSILON);

    return sem;
}

/* 	$Id$	 */


