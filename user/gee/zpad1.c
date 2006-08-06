/* Zero pad.  Surround data by zeros. 1-D */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

void zpad1_lop(bool adj, bool add, int nd, int np, float *data,  float *padd)
{
    int p,d;

    sf_adjnull(adj,add,nd,np,data,padd);

    for (d=0; d < nd; d++) { 
	p = d + (np-nd)/2;
	if (adj) data[d] += padd[p];
	else     padd[p] += data[d];
    }
}
