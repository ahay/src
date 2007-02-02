/* Multiscale helical filter */
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
/*^*/

#include "mshelix.h"

#ifndef _mshelix_h

typedef struct mshelixfilter {
    int nh, ns;
    float* flt;
    int** lag;
    bool** mis;
    sf_filter one;
} *msfilter;
/*^*/

#endif

void onescale(int i, msfilter aa) 
/*< select one scale from multiple scales >*/
{
    aa->one->lag = aa->lag[i];
    if (NULL != aa->mis) aa->one->mis = aa->mis[i];
}

msfilter msallocate(int nh /* filter size */, 
		    int ns /* number of scales */) 
/*< allocate filter >*/
{
    msfilter aa;

    aa = (msfilter) sf_alloc(1,sizeof(*aa));
    aa->nh = nh;
    aa->ns = ns;
    aa->flt = sf_floatalloc(nh);
    aa->lag = sf_intalloc2(nh,ns);
    aa->mis = NULL;
    aa->one = (sf_filter) sf_alloc(1,sizeof(*aa->one));
    aa->one->flt = aa->flt;
    aa->one->nh = nh;

    return aa;
}

void msdeallocate( msfilter aa) 
/*< free allocated storage >*/
{
    free( aa->flt);
    free( aa->lag[0]);
    free( aa->lag);
    if (NULL != aa->mis) {
	free( aa->mis[0]);
	free( aa->mis);
    }
    free( aa->one);
    free( aa);
}

/* 	$Id$	 */
