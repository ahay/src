/* 1-D wave-equation time step using the lowrank approximation */
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

#include "lowrank1.h"

static int nx;
static float *cwave, **wave, *prev;

void lowrank1_init(int n1 /* data size */, 
		   int m1 /* maximum rank */)
/*< initialize >*/
{
    int ix;

    nx = n1;

    cwave = sf_floatalloc(nx);
    prev = sf_floatalloc(nx);
    wave = sf_floatalloc2(nx,m1);

    for (ix = 0; ix < nx; ix++) {
	prev[ix] = 0;
    }

    sf_cosft_init(nx);
}

void lowrank1_close(void)
/*< free allocated storage >*/
{
    free(cwave);
    free(prev);
    free(*wave);
    free(wave);
}

void lowrank1_step(int m1, int m2, 
		   float **lft, 
		   float **mid, 
		   float **rht, 
		   float* curr)
/*< time step >*/
{
    int ix,ik,im;
    float *wavem, old, f;
    
    for (ix = 0; ix < nx; ix++) {
	cwave[ix] = curr[ix];
    }
    sf_cosft_frw(cwave,0,1);
    
    for (im = 0; im < m2; im++) {
	wavem = wave[im];
	for (ik = 0; ik < nx; ik++) {
	    wavem[ik] = cwave[ik]*rht[ik][im];
	}
	sf_cosft_inv(wavem,0,1);
    }
    
    for (ix = 0; ix < nx; ix++) {
	old = f = curr[ix];
	f += f-prev[ix];
	prev[ix] = old;
	
	for (im = 0; im < m2; im++) {
	    for (ik = 0; ik < m1; ik++) {
		f += lft[ik][ix]*mid[im][ik]*wave[im][ix];
	    }
	}
	curr[ix] = f;
    } 
}
