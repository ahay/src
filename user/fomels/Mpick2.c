/* Automatic picking  from semblance-like panels. 

Takes: rect1=1 rect2=1 ...

rectN defines the size of the smoothing stencil in N-th dimension.
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
#include <rsf.h>

#include "dynprog2.h"
#include "divn.h"

int main(int argc, char* argv[])
{
    int dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    int it, niter, nm, n1, n2, n3, n4, i4, i3, i2, i1, i, gate, i0, **ipick;
    float ***scan, ***weight, *pick, *ampl, *pick2, o2, d2, an, asum, a, vel0;
    bool smooth;
    char key[6], *label;
    sf_file scn, pik;

    sf_init(argc,argv);
    scn = sf_input("in");
    pik = sf_output("out");

    if (SF_FLOAT != sf_gettype(scn)) sf_error("Need float input");
    dim = sf_filedims (scn,n);
    if (dim < 3) sf_error("Need at least three dimensions");

    n4 = 1;
    for (i=3; i < dim; i++) {
	n4 *= n[i];
    }

    n1 = n[0];
    n2 = n[1];
    n3 = n[2];
    nm = n1*n3*n4;

    if (!sf_histfloat(scn,"o2",&o2)) o2=0.;
    if (!sf_histfloat(scn,"d2",&d2)) d2=1.;
 
    if (!sf_getfloat("vel0",&vel0)) vel0=o2;
    /* surface velocity */
    i0 = 0.5 + (vel0-o2)/d2;
    if (i0 < 0) i0=0;
    if (i0 >= n2) i0=n2-1;

    sf_putint(pik,"n2",1);
    if (NULL != (label = sf_histstring(scn,"label2"))) 
	sf_putstring(pik,"label",label);
    if (NULL != (label = sf_histstring(scn,"unit2"))) 
	sf_putstring(pik,"unit",label);

    for (i=1; i < dim-1; i++) {
	n[i] = n[i+1];
    }
    dim--;

    for (i=0; i < dim; i++) {
	if (n[i] > 1) {
	    snprintf(key,6,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=1;
	} else {
	    rect[i]=1;
	}
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getfloat("an",&an)) an=1.; 
    /* axes anisotropy */
    if (!sf_getint("gate",&gate)) gate=3; 
    /* picking gate */
    if (!sf_getbool("smooth",&smooth)) smooth=true;
    /* if apply smoothing */

    scan = sf_floatalloc3(n1,n2,n3);
    weight = sf_floatalloc3(n2,n3,n1);
    
    dynprog2_init(n1,n2,n3,gate,gate,an,1.);

    ipick = sf_intalloc2(n1,n3);

    if (smooth) {
	pick = sf_floatalloc(nm);
	pick2 = sf_floatalloc(nm);	
	ampl = sf_floatalloc(nm);

	divn_init(dim,nm,n,rect,niter);
    } else {
	pick = NULL;
	pick2 = sf_floatalloc(n1);
	ampl = NULL;
    }

    for (i4=0; i4 < n4; i4++) {
	sf_warning("record %d of %d",i4+1,n4);
	
	sf_floatread(scan[0][0],n1*n2*n3,scn);

	/* transpose and reverse */
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    weight[i1][i3][i2] = expf(-scan[i3][i2][i1]);
		}
	    }
	}

	dynprog2(i0, weight);
	dynprog2_traj(ipick);

	for (i3=0; i3 < n3; i3++) {
	    if (smooth) {
		for (i1=0; i1 < n1; i1++) {
		    i = i1 + n1*(i3+i4*n3);
		    it = ipick[i3][i1];
		    pick[i] = it;
		    ampl[i] = scan[i1][it][i3];
		}
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    pick2[i1] = o2+ipick[i3][i1]*d2;
		}
		sf_floatwrite(pick2,n1,pik);
	    }
	}

    
	if (smooth) {
	    /* normalize amplitudes */
	    asum = 0.;
	    for (i = 0; i < nm; i++) {
		a = ampl[i];
		asum += a*a;
	    }
	    asum = sqrtf (asum/nm);
	    for(i=0; i < nm; i++) {
		ampl[i] /= asum;
		pick[i] = (o2+pick[i]*d2)*ampl[i];
	    }
    
	    divn(pick,ampl,pick2);

	    sf_floatwrite(pick2,nm,pik);
	} 
    }

    exit(0);
}
