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
    bool segy, desc;
    int i1, i2, n1, n2;
    double *max=NULL, *min=NULL, *mean=NULL;
    int *inp=NULL, *indxmax=NULL, *indxmin=NULL;
    float *inpfloat=NULL;
    double value;
    sf_file head=NULL;
    sf_datatype typehead;

    sf_init (argc,argv);
    head = sf_input("in");

    typehead = sf_gettype(head);
    if (SF_INT != typehead && SF_FLOAT != typehead ) sf_error("Need float or int input headers");
    if (!sf_histint(head,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(head,1);

    fprintf(stdout,"%d headers, %d traces\n",n1,n2);
 
    if (SF_INT == typehead) inp = sf_intalloc(n1);
    else                    inpfloat = sf_floatalloc(n1);  /* input must be float or int */

    min = (double*) sf_alloc(n1,sizeof(double));
    max = (double*) sf_alloc(n1,sizeof(double));

    indxmin = sf_intalloc(n1);
    indxmax = sf_intalloc(n1);

    mean = (double*) sf_alloc(n1,sizeof(double));

    for (i1=0; i1 < n1; i1++) {
	mean[i1] = 0;
    }

    for (i2=0; i2 < n2; i2++) {
	if (SF_INT == typehead) sf_intread  (inp,n1,head);
	else                    sf_floatread(inpfloat,n1,head);
	for (i1=0; i1 < n1; i1++) {
	    if (SF_INT == typehead) value=inp[i1];
	    else                    value=inpfloat[i1];
	    if (i2==0 || (min[i1] > value)) {
		min [i1] = value;
		indxmin[i1] = i2;
	    }
	    if (i2==0 || (max[i1] < value)) {
		max [i1] = value;
		indxmax[i1] = i2;
	    }
	    mean[i1] += value;
	}
    }

    if (!sf_getbool("segy",&segy)) segy=true;
    /* if SEGY headers */

    if (!sf_getbool("desc",&desc)) desc=false;
    /* if describe keys */

    if (segy) {
	segy_init(n1,head);
    } else {
	other_init(n1,head);
    }
    
    /* put headers on table of numbers */
    printf("******************************************************************************* \n");
    printf("     key     \t            min     \t              max    \t          mean\n");
    printf("------------------------------------------------------------------------------- \n");
    for (i1=0; i1 < n1; i1++) {
	if (min[i1] != 0 || max[i1] != 0) {
	    if (SF_INT == typehead) {
		printf("%-8s %4d %14ld @ %d\t%14ld @ %d\t%14g\n",
		       segykeyword(i1),i1,
		       lrint(min[i1]),indxmin[i1],
		       lrint(max[i1]),indxmax[i1],
		       mean[i1]/n2);
	    } else {
		printf("%-8s %4d %14g @ %d\t%14g @ %d\t%14g\n",
		       segykeyword(i1),i1,
		       min[i1],indxmin[i1],
		       max[i1],indxmax[i1],
		       mean[i1]/n2);
	    }
	    if (desc) printf("[%s]\n",segydesc(i1));
	} 
    }
    printf("******************************************************************************* \n");

    exit(0);
}
