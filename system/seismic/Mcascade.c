/* Velocity partitioning for cascaded migrations. */
/*
  Copyright (C) 2007 University of Texas at Austin

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
int main(int argc, char* argv[])
{
    int n1, n2, i1, i2, *ntcut=NULL, ncut, ic, nt;
    float dt, o1, *tcut=NULL, *vint=NULL, *vtmp=NULL, vt;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1) || 0. != o1)
	sf_error("Need o1=0 in input");

    n2 = sf_leftsize(in,1);

    if (!sf_getint("ncut",&ncut)) ncut=1;
    /* number of cuts */

    sf_putint(out,"n2",ncut+1);
    sf_putfloat(out,"d2",1.);
    sf_putfloat(out,"o2",0.);

    sf_shiftdim(in, out, 2);

    tcut = sf_floatalloc(ncut);
    ntcut = sf_intalloc(ncut);

    if (!sf_getints("ntcut",ntcut,ncut)) {
	if (!sf_getfloats("tcut",tcut,ncut)) 
	    sf_error("Need tcut= or ntcut=");
	/* time cuts */
	for (ic=0; ic < ncut; ic++) {
	    ntcut[ic] = tcut[ic]/dt;
	    sf_warning("ntcut[%d]=%d",ic,ntcut[ic]);
	}
    }

    vint = sf_floatalloc(n1);
    vtmp = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(vint,n1,in);

	for (i1=0; i1 < n1; i1++) {
	    vint[i1] *= vint[i1];
	}

	for (ic=0; ic < ncut; ic++) {
	    nt = ntcut[ic];
	    vt = vint[nt];

	    for (i1=0; i1 < nt; i1++) {
		vtmp[i1] = sqrtf(vint[i1]);
		vint[i1] = 0.;
	    }

	    for (i1=nt; i1 < n1; i1++) {
		vtmp[i1] = sqrtf(vt);
		vint[i1] -= vt;
		if (vint[i1] < 0.) sf_error("Negative velocity at %d cut",ic+1);
	    }

	    sf_floatwrite(vtmp,n1,out);
	}

	for (i1=0; i1 < n1; i1++) {
	    vtmp[i1] = sqrtf(vint[i1]);
	}

	sf_floatwrite(vtmp,n1,out);
    }

    exit(0);
}
