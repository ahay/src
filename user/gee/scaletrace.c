/* Trace scaling. */
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

static int n1, n2;
static const float *data;

void scaletrace_init(int n1_in, int n2_in, const float *data_in)
/*< initialize >*/
{
    n1 = n1_in;
    n2 = n2_in;
    data = data_in;
}

void scaletrace_lop(bool adj, bool add, int ns, int nd, 
		    float *scale, float *sdata)
/*< linear operator >*/
{
    int i1,i2,i;
    
    sf_adjnull (adj, add, ns, nd, scale, sdata);

    for (i=i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++, i++) {
	    if( adj) scale[i2] += sdata[i] * data[i];
	    else     sdata[i] += scale[i2] * data[i];
	}
    }
}

