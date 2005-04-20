/* Bandpass filtering */
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

#include <math.h>
#include <stdio.h>

#include "bandpass.h"
#include "butter.h"

static butter blo, bhi;

void bandpass_init (void)
/*< initialize (using a fixed band definition) >*/
{
    const float fhi=0.275314, flo=0.0141445;
    const int nphi=10, nplo=5;

    blo = butter_init(false,flo,nplo);
    bhi = butter_init(true, fhi,nphi);
}

void bandpass_close (void)
/*< free allocated storage >*/
{
    butter_close(blo);
    butter_close(bhi);
}

void bandpass (int n1, float* trace /* [n1] */)
/*< bandpass a trace (in place) >*/
{
    int i1;
    float t;

    butter_apply(blo,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    butter_apply(blo,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    butter_apply(bhi,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
    butter_apply(bhi,n1,trace);
    for (i1=0; i1 < n1/2; i1++) { 
	t=trace[i1];
	trace[i1]=trace[n1-1-i1];
	trace[n1-1-i1]=t;
    }
}

/* 	$Id$	 */
