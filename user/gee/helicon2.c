/* Double convolution with a helix filter */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "helicon2.h"
#include "helicon.h"

static float* tmp;

void helicon2_init(int nd /* data size */,
		   sf_filter bb) 
/*<  Initialized with the filter. >*/
{
    tmp = sf_floatalloc(nd);
    helicon_init(bb);
}

void helicon2_lop( bool adj, bool add, 
		   int nx, int ny, float* xx, float*yy) 
/*< linear operator >*/
{
    sf_chain(helicon_lop,
	     helicon_lop,
	     adj,add,nx,ny,nx,xx,yy,tmp);
}

void helicon2_close(void)
/*< free allocated storage >*/
{
    free(tmp);
}
