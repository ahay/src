/* Make synthetic topography map. */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

#include "random.h"

int main(int argc, char* argv[])
{
    int n1, n2, nm, i1, i2, *known;
    float *hh, *m1, *m2, *s1, *s2;
    sf_triangle t1, t2;
    sf_file syn, mod, mask;
 
    sf_init (argc,argv);
    mod = sf_output("mod");
    syn = sf_output("out");
    mask = sf_output("mask");

    if (!sf_getint("n1",&n1)) n1=100;
    if (!sf_getint("n2",&n2)) n2=100;
    /* data dimensions */

    sf_putint(syn,"n1",n1);
    sf_putint(syn,"n2",n2);
    sf_setformat(syn,"native_float");

    sf_putint(mod,"n1",n1);
    sf_putint(mod,"n2",n2);
    sf_setformat(mod,"native_float");

    sf_putint(mask,"n1",n1);
    sf_putint(mask,"n2",n2);
    sf_setformat(mask,"native_int");

    nm = 35*n2;

    hh = sf_floatalloc(n1);
    m1 = sf_floatalloc(nm);
    m2 = sf_floatalloc(nm);
    s1 = sf_floatalloc(n1);
    s2 = sf_floatalloc(n2);
    known = sf_intalloc(n1);

    random_init (1993);

    for (i1=0; i1 < nm; i1++) {
	m1[i1] = powf(random0()-0.5,7);
	m2[i1] = powf(random0()-0.5,7);
    }

    random_init (1993);

    for (i2=0; i2 < n2; i2++) {
	s2[i2] = random0();
    }
    for (i1=0; i1 < n1; i1++) {
	s1[i1] = random0();
    }

    t1 = sf_triangle_init(2,nm);
    sf_smooth2(t1,0,1,false,false,m1);

    t2 = sf_triangle_init(50,nm);
    sf_smooth2(t2,0,1,false,false,m2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    hh[i1] = 
		m1[i2+1 +2*n2 + i1] +
		m2[i2  +17*n2 - (i1+1)*15] * 10. +
		cosf((i2+1)*.5 - (i1+1)*.75) / 500.;
	}
	sf_floatwrite(hh,n1,mod);
	for (i1=0; i1 < n1; i1++) {
	    if ( i1*i1 + i2*i2 > (n1*n1 + n2*n2)/8 &&
		 s2[i2] >.15 && s1[i1] >.15) {
		known[i1] = 0;
		hh[i1] = 0.;
	    } else {
		known[i1] = 1;
	    }
	}
	sf_floatwrite(hh,n1,syn);
	sf_intwrite(known,n1,mask);
    }

    exit(0);
}
