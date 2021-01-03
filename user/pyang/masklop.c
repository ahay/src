/* mask linear operator for 2D test
 */
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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

#include "masklop.h"

static int n1, n2;
static float *mask;

void mask_init(int n1_, int n2_, float *mask_)
/*< initialize mask operator >*/
{
    n1=n1_;
    n2=n2_;
    mask=mask_;
}

void mask_lop(bool adj, bool add, int nm, int nd, sf_complex *mm, sf_complex *dd)
/*< linear operator >*/
{
    int i1, i2;
    if(nm!=n1*n2) sf_error("datasize mismatch");
    if(nd!=n1*n2) sf_error("datasize mismatch");

    sf_cadjnull(adj, add, nm, nd, mm, dd);

    if(adj){
	for(i2=0; i2<n2; i2++){
	    if(mask[i2]){
		for(i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H	    
		    mm[i1+n1*i2]+=dd[i1+n1*i2];
#else
		    mm[i1+n1*i2]=sf_cadd(mm[i1+n1*i2],dd[i1+n1*i2]);
#endif
		}
	    }
	}
    }else{
	for(i2=0; i2<n2; i2++){
	    if(mask[i2]){
		for(i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H	  
		    dd[i1+n1*i2]+=mm[i1+n1*i2];
#else
		    dd[i1+n1*i2]=sf_cadd(dd[i1+n1*i2],mm[i1+n1*i2]);
#endif
		}
	    }
	}
    }
}
