/* Automatic picking from semblance-like panels (3-D input). */

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
#include <rsf.h>

#include "dynprog.h"

int main(int argc, char* argv[])
{
    int it, niter, nm, n1, n2, n3, i3, i2, i1, i, gate, i0, rect1, rect2, n123, n[2], rect[2];
    float **scan, **weight, *pick, *ampl, *pick2, o2, d2, an, asum, a, ct, vel0;
    char *label;
    sf_file scn, pik;

    sf_init(argc,argv);
    scn = sf_input("in");
    pik = sf_output("out");

    if (SF_FLOAT != sf_gettype(scn)) sf_error("Need float input");

    if (!sf_histint(scn,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(scn,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(scn,"n3",&n3)) sf_error("No n3= in input");

    nm = n1*n3;
    n123 = nm*n2;
    n[0]=n1;
    n[1]=n3;

    if (!sf_histfloat(scn,"o2",&o2)) o2=0.;
    if (!sf_histfloat(scn,"d2",&d2)) d2=1.;
 
    if (!sf_getfloat("vel0",&vel0)) vel0=o2;
    /* surface velocity */
    i0 = 0.5 + (vel0-o2)/d2;
    if (i0 < 0) i0=0;
    if (i0 >= n2) i0=n2-1;

    sf_unshiftdim(scn,pik,2);

    if (NULL != (label = sf_histstring(scn,"label2"))) 
		sf_putstring(pik,"label",label);
    if (NULL != (label = sf_histstring(scn,"unit2"))) 
		sf_putstring(pik,"unit",label);

    if (!sf_getint("rect1",&rect1)) rect1=1;
    /* smoothing radius on the first axis */ 
    if (!sf_getint("rect2",&rect2)) rect2=1;
    /* smoothing radius on the second axis */ 
    rect[0]=rect1;
    rect[1]=rect2;

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getfloat("an",&an)) an=1.; 
    /* axes anisotropy */
    if (!sf_getint("gate",&gate)) gate=3; 
    /* picking gate */

    scan = sf_floatalloc2(n1,n2);
    weight = sf_floatalloc2(n2,n1);

    (void) dynprog_init(n1,n2,gate,an,false);

    pick = sf_floatalloc(nm);
    pick2 = sf_floatalloc(nm);	
    ampl = sf_floatalloc(nm);
    
    sf_divn_init(2,nm,n,rect,niter,true);

    for (i3=0; i3 < n3; i3++) {
	sf_warning("cmp %d of %d;",i3+1,n3);
	sf_floatread(scan[0],n1*n2,scn);

	/* transpose and reverse */
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		weight[i1][i2] = expf(-scan[i2][i1]);
	    }
	}

	dynprog(i0, weight);
	dynprog_traj(pick2);

	for (i1=0; i1 < n1; i1++) {
	    i = i1 + i3*n1;
	    ct = pick2[i1];
	    pick[i] = ct;
	    it = floorf(ct);
	    ct -= it;
	    if (it >= n2-1) {
		ampl[i]=scan[n2-1][i1];
	    } else if (it < 0) {
		ampl[i]=scan[0][i1];
	    } else {
		ampl[i]=scan[it][i1]*(1.-ct)+scan[it+1][i1]*ct;
	    }
	}
    }
    sf_warning(".");
    
    /* normalize amplitudes */
    asum = 0.;
    for (i = 0; i < nm; i++) {
	a = ampl[i];
	asum += a*a;
    }
    asum = sqrtf (asum/nm);
    for(i=0; i < nm; i++) {
	ampl[i] /= asum;
	pick[i] = (o2+pick[i]*d2-vel0)*ampl[i];
    }
    
    sf_divn(pick,ampl,pick2);

    for(i=0; i < nm; i++) {
	pick2[i] += vel0;
    }

    sf_floatwrite(pick2,nm,pik);

    exit(0);
}
