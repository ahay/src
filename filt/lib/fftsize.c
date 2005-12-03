/* Select optimal size for KISS FFT */
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
#include <float.h>

#include "fftsize.h"
#include "_fftsize.h"

int sf_fftr_size(int min)
/*< Find optimal size for real FFT greater or equal min >*/
{
    int n, size;

    for (n=0; n < SF_FFTR_SIZE; n++) {
	size = SF_FFTR[n].size;
	if (size >= min) return size;
    }
    return min;
}

int sf_fftr_size2(int min, int max)
/*< Find optimal size for real FFT between min and max >*/
{
    int n, size, nmin, nmax;
    float cost, mincost;

    nmin=-1;
    for (n=0; n < SF_FFTR_SIZE; n++) {
	size = SF_FFTR[n].size;
	if (size >= min) {
	    nmin = n;
	    break;
	}
    }
    if (nmin < 0) return min;

    nmax=0;
    for (n=SF_FFTR_SIZE-1; n >= 0; n--) {	
	size = SF_FFTR[n].size;
	if (size <= max) {
	    nmax = n;
	    break;
	}
    }

    mincost = FLT_MAX;
    for (n=nmin; n <= nmax; n++) {
	cost = SF_FFTR[n].cost;
	if (cost < mincost) {
	    mincost = cost;
	    size=SF_FFTR[n].size;
	}
    }
    return size;
}
