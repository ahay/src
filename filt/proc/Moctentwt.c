/* Tent-like weight for out-of-core patching. */
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

#include <unistd.h>

#include <rsf.h>

#include "ocmkwallwt.h"
#include "oc.h"

int main(int argc, char* argv[])
{
    int dim, j, *n, *w, *k, *a, *center, n12, w12;
    float *wind;
    FILE *wall;
    char string[256], *temp;
    sf_file wallwt, windwt;

    sf_init (argc, argv);
    wallwt = sf_output("out");
    windwt = sf_output("windwt");

    sf_setformat(wallwt,"native_float");
    sf_setformat(windwt,"native_float");

    if (!sf_getint ("dim",&dim)) dim = 2;
    /* number of dimensions */

    n = sf_intalloc (dim);
    w = sf_intalloc (dim);
    k = sf_intalloc (dim);
    center = sf_intalloc (dim);
    a = sf_intalloc (dim);

    n12 = 1;
    for(j=0; j < dim; j++) {
	sprintf(string,"n%d",j+1);
	if (!sf_getint(string,&n[j])) n[j]=1;
	sf_putint(wallwt,string,n[j]);
	n12 *= n[j];
    }

    if (!sf_getints ("w",w,dim)) sf_error("Need w=");
    /* window size */

    w12 = 1;
    for(j=0; j < dim; j++) {
	if (w[j] > n[j]) w[j] = n[j];
	k[j] = 1.5 * n[j] / (w[j] - a[j] + 1.);
	a[j] = 1;
	center[j] = 1;
	w12 *= w[j];
	sprintf(string,"n%d",j+1);
	sf_putint(windwt,string,w[j]);
    }

    sf_getints ("k",k,dim);
    /* number of windows */

    sf_getints ("a",a,dim);
    /* filter size */
    
    for(j=1; j < dim; j++) {
	if (a[j] > 1) center[j-1] = a[j-1]/2;
    }

    sf_getints ("center",center,dim);

    wind = sf_floatalloc(w12);
    wall = sf_tempfile(&temp,"w+b");

    sf_tent2 (dim, w, wind);

    sf_floatwrite (wind, w12, windwt);

    ocmkwallwt (false, dim, k, n, w, wind, wall);

    oc_dump(n12*sizeof(float),wall,wallwt);

    unlink(temp);

    exit (0);
}

/* 	$Id: Mtentwt.c 691 2004-07-04 19:28:08Z fomels $	 */
