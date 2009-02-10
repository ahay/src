/* Delete-d jackknife estimator's variance. */
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
#include "combinatory.h"

int main(int argc, char* argv[])
{
    int i,ns,n,d;
    float *sj,tm,sp,sm,*stde;

    sf_axis jstde;
    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("n",&n)) sf_error("No n=");
    /* number of total sample */

    if (!sf_getint("d",&d)) sf_error("No d=");
    /* number of d delete samples */

    if (!sf_histint(in,"n1",&ns)) sf_error("No n1=");
    /* number of jackknife samples (combinations) = n!/(d!(n-d)!) */
    if (ns != binomial(n,d)) sf_error("Incompatible n,d and ns");

    /* Output combinations */
    comb_show(n,d);

   /* output file parameters */
    jstde = sf_maxa(1,0,1);
    sf_oaxa(out,jstde,1);
    stde = sf_floatalloc(1);

    /* memory allocations */
    sj = sf_floatalloc(ns);

    /* read input jackknife samples */
    sf_floatread(sj,ns,in);

    /* Expectation of jackknife replications sj of statistical estimator */
    sp = 0.0;
    for (i = 0; i < ns; i++) {
	sp += sj[i];
    }
    sp /= ns;

    /* Delete-d jackknife estimator's standard error */
    sm = 0.0;
    for (i = 0; i < ns; i++) {
	tm = sj[i]-sp;
	sm += tm*tm;
    }
    stde[0] = sqrtf(sm*(n-d)/(d*ns));

    sf_floatwrite(stde,1,out);
    
    exit(0);
}


