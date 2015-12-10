/* Bandpass filtering */

/*
  Copyright (C) 2015 University of Texas at Austin

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
static sf_butter blo, bhi;

void bandpass_init(int n                /* samples in a trace */,
		   float flo, float fhi /* low and high frequencies */,
		   int nplo, int nphi   /* number of poles */)
/*< initialize >*/
{
    n1 = n;
    blo = sf_butter_init(false, flo, nplo);
    bhi = sf_butter_init(true,  fhi, nphi);
}

void bandpass_close(void)
/*< free allocated storage >*/
{
    sf_butter_close(blo);
    sf_butter_close(bhi);
}

void bandpass(float* trace)
/*< apply bandpass filtering >*/
{
    sf_butter_apply (blo, n1, trace); 
    sf_reverse (n1, trace);
    sf_butter_apply (blo, n1, trace); 
    sf_reverse (n1, trace);
    sf_butter_apply (bhi, n1, trace);
    sf_reverse (n1, trace);
    sf_butter_apply (bhi, n1, trace); 
    sf_reverse (n1, trace);
}
