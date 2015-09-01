/* Automatic picking  from semblance-like panels. 

Takes: rect1=1 rect2=1 ...

rectN defines the size of the smoothing stencil in N-th dimension.

Theory in Appendix B of:
S. Fomel, 2009, 
Velocity analysis using AB semblance: Geophysical Prospecting, v. 57, 311-321.
Reproducible version in RSFSRC/book/tccs/avo 
http://ahay.org/RSF/book/tccs/avo/paper_html/

August 2012 program of the month:
http://ahay.org/blog/2012/08/01/program-of-the-month-sfpick/
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

#include "dynprog.h"

int main(int argc, char* argv[])
{
    int dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    int it, niter, nm, n1, n2, n3, i3, i2, i1, i, gate, i0;
    float **scan, **weight, *pick, *ampl, *pick2, o2, d2, an, asum, a, ct, vel0;
    bool smooth;
    char key[6], *label;
    sf_file scn, pik;

    sf_init(argc,argv);
    scn = sf_input("in");
    pik = sf_output("out");

    if (SF_FLOAT != sf_gettype(scn)) sf_error("Need float input");
    dim = sf_filedims (scn,n);
    if (dim < 2) sf_error("Need at least two dimensions");

    n3 = 1;
    for (i=2; i < dim; i++) {
	n3 *= n[i];
    }

    n1 = n[0];
    n2 = n[1];
    nm = n1*n3;

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

    for (i=1; i < dim-1; i++) {
	n[i] = n[i+1];
    }
    dim--;

    for (i=0; i < dim; i++) {
	if (n[i] > 1) {
	    snprintf(key,6,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=1;
	    /*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
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

    scan = sf_floatalloc2(n1,n2);
    weight = sf_floatalloc2(n2,n1);

    (void) dynprog_init(n1,n2,gate,an,false);

    if (smooth) {
	pick = sf_floatalloc(nm);
	pick2 = sf_floatalloc(nm);	
	ampl = sf_floatalloc(nm);

	sf_divn_init(dim,nm,n,rect,niter,true);
    } else {
	pick = NULL;
	pick2 = sf_floatalloc(n1);
	ampl = NULL;
    }

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

	if (smooth) {
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
	} else {
	    for (i1=0; i1 < n1; i1++) {
		pick2[i1] = o2+pick2[i1]*d2;
	    }
	    sf_floatwrite(pick2,n1,pik);
	}
    }
    sf_warning(".");
    
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
	    pick[i] = (o2+pick[i]*d2-vel0)*ampl[i];
	}
    
	sf_divn(pick,ampl,pick2);

	for(i=0; i < nm; i++) {
	    pick2[i] += vel0;
	}

	sf_floatwrite(pick2,nm,pik);
    } 

    exit(0);
}
