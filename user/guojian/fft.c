/* FFT routines. */
/*
  Copyright (C) 2006 The Board of Trustees of Stanford University
  
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

static int n1;
static float wt;
static kiss_fft_cfg cfg1=NULL, icfg1=NULL;

#include "fft.h"

void kissfft(void *fz, int n, int isign)
/*< 1-D complex FFT >*/
{ 
  int i; 
  kiss_fft_cpx *cz = (kiss_fft_cpx *) fz;
  sf_complex *cc = (sf_complex *) fz;
  
  if (isign > 0) {
      for (i=0; i<n1; i++) { /* FFT centering */
	  cc[i] = (i%2? 1:-1) * cc[i];
      }      
      kiss_fft_stride(cfg1,cz,cz,1);
  } else {
      kiss_fft_stride(icfg1,cz,cz,1);
      for (i=0; i<n1; i++) {
	  cc[i] = (i%2? wt:-wt) * cc[i];
      }
  }
}

int kissfftn(int nmin)
/*< next size >*/
{
    n1 = kiss_fft_next_fast_size(nmin);
    cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
    icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);
    wt =  1.0/n1;
    return n1;
}
