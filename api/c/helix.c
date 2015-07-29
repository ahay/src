/* Helical filter definition and allocation. */
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
#include <stdlib.h>

#include "helix.h"
#include "alloc.h"
#include "error.h"

#include "_bool.h"
/*^*/

#ifndef _sf_helix_h

typedef struct sf_helixfilter {
    int     nh;
    float* flt;
    int*   lag;
    bool*  mis;
    float   h0;
} *sf_filter;
/*^*/

#endif

/*------------------------------------------------------------*/
sf_filter sf_allocatehelix( int nh) 
/*< allocation >*/
{
    sf_filter aa;

    aa = (sf_filter) sf_alloc(1,sizeof(*aa));
    aa->nh = nh;
    if (nh > 0) {
	aa->flt = sf_floatalloc(nh);
	aa->lag = sf_intalloc  (nh);
    } else {
	aa->flt = NULL;
	aa->lag = NULL;
    }
    aa->mis = NULL;
    
    return aa;
}

/*------------------------------------------------------------*/
void sf_deallocatehelix( sf_filter aa) 
/*< deallocation >*/
{
    free( aa->flt);
    free( aa->lag);
    if (NULL != aa->mis) free( aa->mis);
    free( aa);
}

/*------------------------------------------------------------*/
void sf_displayhelix( sf_filter aa)
/*< display filter >*/
{
    int ih;

    sf_warning("------------------------------------------------------------");
    sf_warning("h0=%g",aa->h0);
    for (ih=0; ih < aa->nh; ih++) {
	sf_warning("%4d %4d %4.8g",ih,aa->lag[ih],aa->flt[ih]);
    }
    sf_warning("------------------------------------------------------------");
}
