/* Singular value decomposition (SVD)
   Compute [U,O,V]=SVD(A), A=UOV, a little bit different from A=UoV' */
/*
  Copyright (C) 2013 University of Texas at Austin
   
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
/* Author: Yangkang Chen */

#include <rsf.h>

#include "svd.h"

int main(int argc, char *argv[])
{
    int i, j, n1, n2, n3, ka; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    sf_file in, out, outu=NULL, outv=NULL;
    const double eps = 1.0e-5;
    float *u, *v, *a, *d;
    bool ifoutu, ifoutv;
    sf_init(argc,argv);

    in =  sf_input("in");    
    /* in is the input matrix A, with n2 rows and n1 columns */
    out = sf_output("out");
    /* out is diagonal matrix o, with n2 rows and n1 columns */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    ka = SF_MAX(n1,n2)+1;

    sf_putint(out, "n2", 1);

    ifoutu = (NULL!=sf_getstring("left"));
    ifoutv = (NULL!=sf_getstring("right"));

    if(ifoutu)
    {
	outu= sf_output("left");
	/* outu is square matrix u,  with n2 rows and n2 columns */
	sf_putint(outu, "n1", n2);  /* outu is u */
	sf_putint(outu, "n2", n2);
	sf_putint(outu, "n3", n3);}

    if(ifoutv)
    {
	outv= sf_output("right");
	/* outv is square matrix v,  with n1 rows and n1 columns */
	sf_putint(outv, "n1", n1);  /* outv is v */
	sf_putint(outv, "n2", n1);
	sf_putint(outv, "n3", n3);}

    svdinit(n2, n1, ka, eps);
    a = sf_floatalloc(n1*n2);
    u = sf_floatalloc(n2*n2);
    v = sf_floatalloc(n1*n1);
    d = sf_floatalloc(n1);

    for(i=0;i<n3;i++)  {
	sf_floatread(a, n1*n2, in);
        svduav(a, u, v); /* implement svd */
	for (j=0; j < n1; j++) {
	    d[j] = a[j*n1+j];
	}
	sf_floatwrite(d, n1, out);   /* output singular values */ 
	if(ifoutu)sf_floatwrite(u, n2*n2, outu);  /* output u */ 
        if(ifoutv)sf_floatwrite(v, n1*n1, outv);  /* output v */
    }

    /* release memory in main function and subroutines */
    svdclose();
    free(a);
    free(u);
    free(v);
    exit(0);
}

