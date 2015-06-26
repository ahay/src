/* Simple 2-D synthetics with crossing plane waves.
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

#include <math.h>

#include <rsf.h>

#include "random.h"

int main(int argc, char* argv[])
{
    int i1, i2, i3, n1, n2, n3, n, ns, p, t1, t2;
    float **pp, *s1, *s2;
    bool second;
    sf_triangle tr1, tr2=NULL;
    sf_file mod;

    sf_init (argc,argv);
    mod = sf_output("out");
    sf_setformat(mod,"native_float");

    if (!sf_getint("n1",&n1)) n1=100;
    if (!sf_getint("n2",&n2)) n2=14;
    if (!sf_getint("n3",&n3)) n3=1;
    /* dimensions */

    if (!sf_getbool("second",&second)) second=true;
    /* if n, only one plane wave is modeled */

    sf_putint(mod,"n1",n1);
    sf_putint(mod,"n2",n2);
    sf_putint(mod,"n3",n3);

    pp = sf_floatalloc2(n1,n2);

    if (!sf_getint("n",&n)) n=3;
    if (!sf_getint("p",&p)) p=3;
    ns = n1*n;
    
    s1 = sf_floatalloc(ns);
    s2 = sf_floatalloc(ns);

    random_init (1993);
    for (i1=0; i1 < ns; i1++) {
		s1[i1] = powf(random0()-0.5,(float) p);
		s2[i1] = powf(random0()-0.5,(float) p);
    }

    if (!sf_getint("t1",&t1)) t1=4;
    /* triangle smoother for first wave */
    tr1 = sf_triangle_init (t1, ns,false);
    sf_smooth2(tr1,0,1,false,s1);

    if (second) {
	if (!sf_getint("t2",&t2)) t2=4;
	/* triangle smoother for second wave */
	tr2 = sf_triangle_init (t2, ns,false);
	sf_smooth2(tr2,0,1,false,s2);
    }	

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    pp[i2][i1]  = s1[i1 + n/2 * n1 + i2];
	    if (second) pp[i2][i1] += s2[i1 + n/2 * n1 - i2*2];
	}
    }

    for (i3=0; i3 < n3; i3++) {
	sf_floatwrite (pp[0],n1*n2,mod);
    }
    
    exit(0);
}

/* 	$Id$	 */
