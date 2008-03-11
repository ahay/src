/* Double causal integration */
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

#include "doubint.h"

static int n;
static float *tt;

void doubint_init(int n1)
/*< initialize >*/
{
    n = n1;
    tt = sf_floatalloc(n);
}

void doubint_close(void)
/*< free allocated storage >*/
{
    free(tt);
}

void doubint_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    if (nx != n || ny != n) sf_error("%s: wrong size",__FILE__);

    sf_adjnull (adj, add, nx, ny, xx, yy);
    sf_chain (sf_causint_lop,
	      sf_causint_lop,
	      adj, true, n, n, n, xx, yy, tt);
}

/* 	$Id: causint.c 2119 2006-08-06 02:03:25Z sfomel $	 */
