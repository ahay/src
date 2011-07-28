/* Monte Carlo integration of cos(2t).P(x1,x2).P(y1,y2) */
/* P is correlated Gaussian joint distribution and t angle between vectors (x1,y1) and (x2,y2) */
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

#include <math.h>
#include <rsf.h>

#include "bmcgauss.h"

int main(int argc, char* argv[])
{
    int i,k;
    int n,nr,iseed;

    float m1,s1;           /* mean and standard deviation */
    float dr,orig,r;         /* correlation coefficient */
    float x1,y1,x2,y2;     /* Gaussian distributed correlated random deviates */
    float cs12,irm1,irm2;

    float *mc;

    sf_file in, out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* correlation coefficient grid */
    if (!sf_histint(in,"n1",&nr)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dr)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&orig)) sf_error("No o1= in input");

    if (!sf_getint("n",&n)) n=100;
    /* number of random deviates pairs */

    if (!sf_getfloat("m1",&m1)) m1=0.0;
    /* mean for deviates */

    if (!sf_getfloat("s1",&s1)) s1=1.0;
    /* standard deviation for deviates */

    if (!sf_getint("iseed",&iseed)) iseed=-33;
    /* random generator seed */
    if (iseed >= 0) sf_error("Need strictly negative iseed");

    /* memory allocations */
    mc = sf_floatalloc(nr);

    for (k = 0; k < nr; k++) {

	r = orig + k*dr;

	mc[k] = 0.0;

	x1 = 0.; x2 = 0.;
	y1 = 0.; y2 = 0.;

        /* Generates n pairs of Gaussian distributed correlated deviates */
	for (i=0; i<n; i++) {
   
            /* Correlated Gaussian deviates */
	    gauss_joint(&iseed,m1,m1,s1,s1,r,&x1,&x2);
	    gauss_joint(&iseed,m1,m1,s1,s1,r,&y1,&y2);

            /* Estimates <cos(2t)> = <2*cos^2(t)-1> */
	    irm1 = 1.0/sqrtf(x1*x1 + y1*y1);
	    irm2 = 1.0/sqrtf(x2*x2 + y2*y2);

	    cs12 = (x1*x2 + y1*y2)*irm1*irm2;
	    mc[k] += (2.*cs12*cs12 - 1.0)/n;

	}
    }

    /* output */
    sf_floatwrite(mc,nr,out);

    exit(0);
}
