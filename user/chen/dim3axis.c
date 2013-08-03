/* 3 dims separated by axis */

/*
  Copyright (C) 2013 Zhonghuan Chen,
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

int dim3axis(sf_file in, int axis, int *n)
/*< 3 dims >*/
{
	int nn[SF_MAX_DIM], dim, i;

	dim  = sf_filedims(in, nn);
	if(dim < axis) return dim;

	n[0] = 1; n[1] = nn[axis-1]; 
	for(i=0; i<axis-1; i++) n[0] *=nn [i];

	n[2] = sf_leftsize(in, axis);

	return dim;	
}

