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

static int opt_size(int min, int tabsize, const sf_Table *table)
{
    int n, size;
    
    for (n=0; n < tabsize; n++) {
	size = table[n].size;
	if (size >= min) {
	    min=size;
	    break;
	}
    }
    return min;
}

int sf_fftr_size(int min)
/*< Find optimal size for real FFT greater or equal min >*/
{
    min = opt_size(min,SF_FFTR_SIZE,SF_FFTR);
    if (min%2) min++; /* make it even */
    return min;
}


int sf_fft_size(int min)
/*< Find optimal size for complex FFT greater or equal min >*/
{
    min = opt_size(min,SF_FFT_SIZE,SF_FFT);
    return min;
}

static int opt_size2(int min, int max, int tabsize, const sf_Table *table)
{
    int n, size, nmin, nmax;
    float cost, mincost;

    size=min;

    nmin=-1;
    for (n=0; n < tabsize; n++) {
	size = table[n].size;
	if (size >= min) {
	    nmin = n;
	    break;
	}
    }
    if (nmin < 0) return min;

    nmax=0;
    for (n=tabsize-1; n >= 0; n--) {	
	size = table[n].size;
	if (size <= max) {
	    nmax = n;
	    break;
	}
    }

    mincost = FLT_MAX;
    for (n=nmin; n <= nmax; n++) {
	cost = table[n].cost;
	if (cost < mincost) {
	    mincost = cost;
	    size=table[n].size;
	}
    }

    return size;
}

int sf_fftr_size2(int min, int max)
/*< Find optimal size for real FFT between min and max >*/
{
    min = opt_size2(min,max,SF_FFTR_SIZE,SF_FFTR);
    if (min%2) min++; /* make it even */
    return min;
}

int sf_fft_size2(int min, int max)
/*< Find optimal size for complex FFT between min and max >*/
{
    min = opt_size2(min,max,SF_FFT_SIZE,SF_FFT);
    return min;
}
