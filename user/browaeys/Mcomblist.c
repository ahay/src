/* Create masks to remove combinations of k elements out of n */
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

#include <float.h>
#include <rsf.h>

#include "comblist.h"

int main(int argc, char* argv[])
{
    int n,k,nc,i,j;

    int done;       /* boolean for more combinations to compute */
    int *a;         /* list of elements in the current combination (not needed at startup) */
    int **mask;

    sf_axis areplic;
    sf_file in,out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    sf_settype(out,SF_INT);

    if (!sf_getint("k",&k)) sf_error("Need k=");
    /* combination of k elements */

    /* input file */
    if (!sf_histint(in,"n1",&n)) sf_error("No n1=");

    nc = binomial(n,k);
    sf_warning("Number of combinations is %3d",nc);

    /* output file parameters */
    areplic = sf_maxa(nc,0,1);
    sf_oaxa(out,areplic,2);

    sf_putstring (out,"label2", "replication");

    /* memory allocations */
    a = sf_intalloc(k);
    mask = sf_intalloc2(n,nc);

    done = 1;
    j = 0;

    while (1) {

        /* Combination of k elements out of n */
	comb_next(n,k,a,&done);

	if (done) break;
        /* done = 1 if more combinations to compute */
        /* done = 0 when the list is exhausted. */

	for (i = 0; i < k; i++) fprintf(stderr," %3d",a[i]);
	fprintf(stderr," \n");

	for (i = 0; i < n; i++) mask[j][i] = 1;
	for (i = 0; i < k; i++) mask[j][a[i]-1] = 0;

	j++;

    }

    /* output */ 
    sf_warning("Number of combinations is %3d",nc);
    sf_intwrite(mask[0],n*nc,out);

    exit(0);
}

