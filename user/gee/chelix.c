/* Complex helical filter definition and allocation. */
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

#include "chelix.h"


#ifndef _chelix_h

typedef struct chelixfilter {
    int     nh;
    sf_complex* flt;
    int*   lag;
    bool*  mis;
} *cfilter;
/*^*/

#endif

/*------------------------------------------------------------*/
cfilter allocatechelix( int nh) 
/*< allocation >*/
{
    cfilter aa;

    aa = (cfilter) sf_alloc(1,sizeof(*aa));
    aa->nh = nh;
    aa->flt = sf_complexalloc(nh);
    aa->lag = sf_intalloc  (nh);
    aa->mis = NULL;
    
    return aa;
}

/*------------------------------------------------------------*/
void deallocatechelix( cfilter aa) 
/*< deallocation >*/
{
    free( aa->flt);
    free( aa->lag);
    if (NULL != aa->mis) free( aa->mis);
    free( aa);
}

/*------------------------------------------------------------*/
void displaychelix( cfilter aa)
/*< display filter >*/
{
    int ih;

    sf_warning("------------------------------------------------------------");
    for (ih=0; ih < aa->nh; ih++) {
	sf_warning("%4d %4d %4.8g %4.8g",ih,aa->lag[ih],crealf(aa->flt[ih]),cimagf(aa->flt[ih]));
    }
    sf_warning("------------------------------------------------------------");
}

void helimakelag(cfilter aa,int nx, int ny)
/*< make lags for a filter >*/
{
    int nh,ih;

    nh=aa->nh;

    if (ny == 1) { /* 2-D case */
	for (ih=0; ih < nh; ih++) {
	    aa->lag[ih] = ih;
	}
    } else {       /* 3-D case */
	for (ih=0; ih < (nh+1)/2; ih++) {
	    aa->lag[ih] = ih; 
	}
	for (; ih < nh; ih++) {
	    aa->lag[ih] = nx-nh+ih;
	}
    }
}
