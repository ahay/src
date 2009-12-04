/* Low-rank matrix approximation (Lexing Ying) */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#ifndef _lowrank_h

typedef void (*samplefunc)(void);
/*^*/

#endif

float **lowrank(int m, int n,    /* matrix size */
		samplefunc func, /* function to sample rows and columns */
		float eps        /* prescribed accuracy */, 
		int npk          /* number of random rows and columns sampled */, 
		int *idx1        /* chosen columns (output) */, 
		int *idx2        /* chosen rows */)
/*< find low-rank matrix approximation >*/
{
    int ir, nr, *rs;
    float **mid=NULL;
    
    nr = npk;

    init_genrand(2009);

    rs = sf_intalloc(nr);

    for (ir=0; ir < nr; ir++) {
	rs[ir] = ceilf((0.5 + genrand_real1())*(m-nr)+ir);
    }

/*    rs = sort(ceil(rand(1,nr)*(m-nr))) + [1:nr]; */

    return mid;
}

