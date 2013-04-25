/* Integer header attributes. 

Only nonzero values are reported.
*/
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

#include <stdio.h>
#include <string.h>
#include <rsf.h>
#include "segy.h"

int main(int argc, char* argv[])
{
    int i1, i2, n1, n2;
    int *max=NULL, *min=NULL, *inp=NULL, *imax=NULL, *imin=NULL;
    double *mean=NULL;
    char pad[] = "                    ", out[21];
    sf_file head=NULL;

    sf_init (argc,argv);
    head = sf_input("in");

    if (SF_INT != sf_gettype(head)) sf_error("Need int input");
    if (!sf_histint(head,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(head,1);

    fprintf(stdout,"******************************************* \n");
    fprintf(stdout,"%d headers, %d traces\n",n1,n2);
 
    inp = sf_intalloc(n1);

    min = sf_intalloc(n1);
    max = sf_intalloc(n1);

    imin = sf_intalloc(n1);
    imax = sf_intalloc(n1);

    mean = (double*) sf_alloc(n1,sizeof(double));

    sf_intread(inp,n1,head);
    for (i1=0; i1 < n1; i1++) {
	min [i1] = inp[i1];
	max [i1] = inp[i1];
	
	imin[i1] = 0;
	imax[i1] = 0;

	mean[i1] = inp[i1];
    }

    for (i2=1; i2 < n2; i2++) {
	sf_intread(inp,n1,head);
	for (i1=0; i1 < n1; i1++) {
	    if (min[i1] > inp[i1]) {
		min [i1] = inp[i1];
		imin[i1] = i2;
	    }
	    if (max[i1] < inp[i1]) {
		max [i1] = inp[i1];
		imax[i1] = i2;
	    }
	    mean[i1] += inp[i1];
	}
    }

    segy_init(n1,head);
    
    /* put headers on table of numbers */
    fprintf(stdout,"\n");
    fprintf(stdout,"indx     key        indx         min       ");
    fprintf(stdout,"indx         max           mean\n");
    fprintf(stdout," key    name         min       value       ");
    fprintf(stdout," max       value          value\n");
    for (i1=0; i1 < n1; i1++) {
	if (min[i1] != 0 || max[i1] != 0) {
	  fprintf(stdout,"%4d %8s %10d %11d %10d %11d %14e\n",
		  i1,segykeyword(i1),
		  imin[i1],min[i1],
		  imax[i1],max[i1],
		  mean[i1]/n2);
	}
    }

    fprintf(stdout,"******************************************* \n");

    exit(0);
}
