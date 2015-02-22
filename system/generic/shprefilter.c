/* Prefilter for shifted linear interpolation */
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
#include <rsf.h>

#include "shprefilter.h"

static float *tmp1, *tmp2; /* temporary storage */
static float shifted[1] = {0.21/0.79};
static float a0=1.26582278; /* normalization */

void shprefilter (int nt, float* dat /* in - data, out - coefficients */)
/*< Convert 1-D data to shifted-linear coefficients >*/
{
    int i;

    /* Initialize filter for converting for c_n*/
    tmp1 = sf_floatalloc (nt);
    tmp2 = sf_floatalloc (nt);
    sf_recfilt_init (nt, 1, shifted);
    
    /* prefilter in n1 direction*/
    for (i = 0; i < nt; i++) {
	tmp1[i] = dat[i];
    }

    sf_recfilt_lop (false, false,nt,nt,tmp1,tmp2);

    for (i = 0; i < nt; i++) {
	dat[i] = tmp2[i]*a0;
    }
	    
    free (tmp1);
    free (tmp2);
    sf_recfilt_close();
}

void shprefilter2d (int nt1, int nt2, float** dat /* in - data, out - coefficients */)
/*< Convert 2-D data to shifted-linear coefficients >*/
{
    int i, j, nt;

    /* Initialize filter for converting for c_n*/
    nt = nt1;
    tmp1 = sf_floatalloc (nt);
    tmp2 = sf_floatalloc (nt);
    sf_recfilt_init (nt, 1, shifted);
    
    /* prefilter in n1 direction*/
    for (j = 0; j < nt; j++) {
	    for (i = 0; i < nt; i++) {
		tmp1[i] = dat[j][i];
	    }

    sf_recfilt_lop (false, false,nt,nt,tmp1,tmp2);

	    for (i = 0; i < nt; i++) {
		dat[j][i] = tmp2[i]*a0;
	    }
    }
    free (tmp1);
    free (tmp2);
    sf_recfilt_close();
    
    nt = nt2;
    tmp1 = sf_floatalloc (nt);
    tmp2 = sf_floatalloc (nt);
    sf_recfilt_init (nt, 1, shifted);
    
    /* prefilter in n2 direction*/
    for (j = 0; j < nt; j++) {
	    for (i = 0; i < nt; i++) {
		tmp1[i] = dat[i][j];
	    }
    sf_recfilt_lop (false, false,nt,nt,tmp1,tmp2);

	    for (i = 0; i < nt; i++) {
		dat[i][j] = tmp2[i]*a0;
	    }
    }

}

      
/* 	$Id: prefilter.c 7107 2011-04-10 02:04:14Z ivlad $	 */
