/* Simple synthetic marine data example. */
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

int main(int argc, char* argv[]) {
    int nt, nh, ny, nz, it, iy, ih, is, iz;
    float *data, **refl, *depth, layer;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");

    if (!sf_getint("nt",&nt)) nt=250;
    /* time samples */
    if (!sf_getint("nh",&nh)) nh=48;
    /* offset samples */
    if (!sf_getint("ny",&ny)) ny=10;
    /* midpoint samples */
    if (!sf_getint("nz",&nz)) nz=25;
    /* depth samples */

    sf_putint(out,"n1",nt);
    sf_putint(out,"n2",nh);
    sf_putint(out,"n3",ny);

    data = sf_floatalloc(nt);
    refl = sf_floatalloc2(nz,ny);
    depth = sf_floatalloc(nz);

    random_init (1992);
    for (iz=0; iz < nz; iz++) {  
	depth[iz] = random0()*nt; /* reflector depth */
	layer =  2.*random0()-1.; /* reflector strength */
	for (iy=0; iy < ny; iy++) {
            /* texture on layer */
	    refl[iy][iz] = (random0()+1.)*layer;   
	}
    }

    for (is=0; is < ny; is++) {     /* shots */
	for (ih=0; ih < nh; ih++) { /* down cable h = (g-s)/2 */
	    for (it=0; it < nt; it++) {
		data[it] = 0.;
	    }
	    for (iz=0; iz < nz; iz++) { /* Add hyperbola */
		iy = (ny-1-is)+(ih-1);  /* y = midpoint */
		iy = (iy+1)%ny;         /* periodic midpoint */
		it = hypotf(depth[iz],5.*ih);   
		if (it < nt) data[it] += refl[iy][iz];
	    }
	    sf_floatwrite(data,nt,out);
	}
    }

    exit(0);
}

/* 	$Id$	 */
