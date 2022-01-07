/* Threshold float/complex inputs given a constant/varying threshold level (Interface).

Methods available:
- soft
- hard
- non-negative Garrote (nng)

Written by: Gilles Hennenfent & Colin Russell, UBC
Created: February 2006
*/
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

static void thrsample(float *in, float *out, bool complex_data, 
		      const char* mode, float thr);

static bool complex_data;
static int ntr, nt;
static float thr;
static const char* mode;

void thr_init(bool fcomplex_data /* if input is complex */,
              int fntr /* number of traces */,
              int fnt /* number of time samples */,
	      float fthr /* thresholding value */,
	      const char* fmode /*thresholding mode */)	      
/*< initialize >*/
{

    complex_data = fcomplex_data;
    complex_data = false;
    ntr = fntr;
    nt = fnt;
    thr = fthr;
    mode = fmode;

}

void thr_lop(float *x, float *y)
/*< apply thresholding >*/
{

    float thrvalue;

    thrvalue = thr;
    
    thrsample(x,y,complex_data,mode,thrvalue);

}

static void thrsample(float *in, float *out, bool complex_data, 
		      const char* mode, float thr)
{
    float isample=0.0;
    float osample=0.0;
    int i;
    
    for (i=0; i<nt*ntr; i++){

	isample = in[i];
	
    	if (fabsf(isample)>thr){
		if (strcmp(mode,"soft") == 0){
			if (isample>thr)
				osample = isample-thr;
			else
				osample = isample+thr;
		}
		if (strcmp(mode,"hard") == 0)
			osample = isample;
		if (strcmp(mode,"nng") == 0)
		osample = isample-(thr*thr/ isample);
    	} else osample=0.0;

	out[i] = osample;

    }	
	
}


