/* IWI interface for dip estimation */
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

#include "dip3.h"
#include "iwidip.h"

static int n[3], n123, prect[3];
static float d[3];
static int porder, pniter, pliter;
static float *pl, *pr;

void iwidip_init(int n1, int n2, int nh,
		 float d1, float d2,
		 int rect1, int rect2, int rect3,
		 int order, int niter, int liter)
/*< initialization >*/
{
    n[0] = n1; n[1] = n2; n[2] = 2*nh+1;
    n123 = n[0]*n[1]*n[2];

    d[0] = d1; d[1] = d2; d[2] = d2;

    prect[0] = rect1; prect[1] = rect2; prect[2] = rect3;

    porder = order; pniter = niter; pliter = liter;

    /* allocate memory */
    pl = sf_floatalloc(n123);
    pr = sf_floatalloc(n123);
}

void iwidip_free()
/*< free >*/
{
    free(pl); free(pr);
}

void iwidip_both(float *image, float *dip)
/*< estimate dip from both directions >*/
{
    int i;

    dip3_init(n[0],n[1],n[2], prect,pliter,false);

    /* left->right */    
    for(i=0; i < n123; i++) {
	pl[i] = 0.;
    }    
    dip3(false, 2,pniter,porder,1,false, image,pl, NULL,-FLT_MAX,+FLT_MAX);
    
    /* right->left */
    for(i=0; i < n123; i++) {
	pr[i] = 0.;
    }    
    dip3(true,  2,pniter,porder,1,false, image,pr, NULL,-FLT_MAX,+FLT_MAX);

    /* average */
    for(i=0; i < n123; i++) {
	dip[i] = 0.5*pl[i]-0.5*pr[i];
    }

    dip3_close();
}
