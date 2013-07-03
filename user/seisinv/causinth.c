/* average alonge offset axis */
/*
  Copyright (C) 2012 China University of Petroleum
  
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

#include "causinth.h"

#include <rsf.h>
/*^*/

static int nt,nm,nh;

void sf_causinth_init(int n1, int n2, int n3)
/*< initialization >*/
{
   nt=n1;
   nm=n2;
   nh=n3;
}

void sf_causinth_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int i1,i2,i3,i; 
    float t;

    sf_adjnull (adj, add, nx, ny, xx, yy);

    for (i2=0;i2<nm;i2++) {
        for (i1=0;i1<nt;i1++) {
               t=0.;
               if (adj) {
                   for (i3=nh-1,i=1;i3>=0;i3--,i++) {
                       t += yy[i3*nm*nt+i2*nt+i1];
                       xx[i3*nm*nt+i2*nt+i1] += t/i;
                    }
                } else {
                   for (i3=0,i=1;i3<=nh-1;i3++,i++) {
                       t += xx[i3*nm*nt+i2*nt+i1];
                       yy[i3*nm*nt+i2*nt+i1] += t/i;
                    }
                }
         }
    }
}

/* 	$Id: causint.c 3588 2008-05-13 13:18:01Z sfomel $	 */
