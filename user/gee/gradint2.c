/* 2-D interpolation */
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

static float *tmp;

void  gradint2_init (float** coord          /* coordinates [nd][2] */, 
		     float o1, float o2, 
		     float d1, float d2,
		     int   n1, int   n2     /* axes */, 
		     sf_interpolator interp /* interpolation function */, 
		     int nf                 /* interpolator length */, 
		     int nd                 /* number of data points */)
/*< initialize >*/
{
    tmp = sf_floatalloc(nd);
    sf_int2_init(coord,o1,o2,d1,d2,n1,n2,interp,nf,nd);
}

void gradint2_lop (bool adj, bool add, int nm, int ny, float* x, float* ord)
/*< linear operator >*/
{ 
    sf_adjnull (adj,add,nm,ny,x,ord);
    sf_chain(sf_igrad1_lop,sf_int2_lop,adj,true,nm,ny,ny,x,ord,tmp); 
}

void gradint2_close (void)
/*< free allocated storage >*/
{
    free(tmp);
    sf_int2_close();
}
