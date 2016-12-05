/* Automatic picking from 3D semblance-like panels plus additional axis. */
/*
  Copyright (C) 2016 Jilin University
  
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

#include "dynprog3.h"

int main(int argc, char* argv[])
{
    int dim, n[SF_MAX_DIM], rect;
    int it0, it1, niter, n1, n2, i2, n3, i3, i1, gate1, gate2, k2, k3, i5, n5;
    float ***scan, ***weight, **pick, *ampl, **pick2;
    float o2, d2, o3, d3, an1, an2, asum, a, ct0, ct1, vel2, vel3;
    bool smooth;
    sf_file scn, pik;

    sf_init(argc,argv);
    scn = sf_input("in");
    pik = sf_output("out");

    if (SF_FLOAT != sf_gettype(scn)) sf_error("Need float input");
    dim = sf_filedims (scn,n);
    if (dim != 3) sf_error("Need three dimensions");

    n1 = n[0];
    n2 = n[1];
    n3 = n[2];

    n5 = sf_leftsize(scn,3);

    if (!sf_histfloat(scn,"o2",&o2)) o2=0.;
    if (!sf_histfloat(scn,"d2",&d2)) d2=1.;

    if (!sf_histfloat(scn,"o3",&o3)) o3=0.;
    if (!sf_histfloat(scn,"d3",&d3)) d3=1.;
 
    if (!sf_getfloat("vel1",&vel2)) vel2=o2;
    if (!sf_getfloat("vel2",&vel3)) vel3=o3;
    /* surface velocity */

    k2 = 0.5 + (vel2-o2)/d2;
    if (k2 < 0) k2=0;
    if (k2 >= n2) k2=n2-1;

    k3 = 0.5 + (vel3-o3)/d2;
    if (k3 < 0) k3=0;
    if (k3 >= n3) k3=n3-1;

    sf_putint(pik,"n2",1);
    sf_putint(pik,"n3",1);
    sf_putint(pik,"n4",2);
    sf_putint(pik,"n5",n5);

    
    if (!sf_getint("rect1",&rect)) rect=1;
    /* smoothing radius */

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getfloat("an1",&an1)) an1=1.;
    if (!sf_getfloat("an2",&an2)) an2=1.;
    /* axes anisotropy */
    if (!sf_getint("gate1",&gate1)) gate1=3;
    if (!sf_getint("gate2",&gate2)) gate2=3;
    /* picking gate */
    if (!sf_getbool("smooth",&smooth)) smooth=true;
    /* if apply smoothing */

    scan = sf_floatalloc3(n1,n2,n3);
    weight = sf_floatalloc3(n2,n3,n1);

    (void) dynprog3_init(n1,n2,n3,gate1,gate2,an1,an2,false);

    if (smooth) {
	pick = sf_floatalloc2(n1,2);
	pick2 = sf_floatalloc2(n1,2);	
	ampl = sf_floatalloc(n1);

	sf_divn_init(1,n1,&n1,&rect,niter,true);
    } else {
	pick = NULL;
	pick2 = sf_floatalloc2(n1,2);
	ampl = NULL;
    }

    for(i5=0; i5 < n5; i5++) {

	sf_floatread(scan[0][0],n1*n2*n3,scn);
	
	/* transpose and reverse */
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    weight[i1][i3][i2] = expf(-scan[i3][i2][i1]);
		}
	    }
	}
	
	dynprog3(k2, k3, weight);
	dynprog3_traj(pick2);
	
	if (smooth) {
	    for (i1=0; i1 < n1; i1++) {
		ct0 = pick2[0][i1];
		it0 = floorf(ct0);
		ct0 -= it0;
		
		if (it0 >= n2-1) {
		    it0 = n2-2;
		    ct0 = 0.;
		} else if (it0 < 0) {
		    it0 = 0;
		    ct0 = 0.;
		}
		
		ct1 = pick2[1][i1];
		it1 = floorf(ct1);
		ct1 -= it1;
		
		if (it1 >= n3-1) {
		    it1 = n3-2;
		    ct1 = 0.;
		} else if (it1 < 0) {
		    it1 = 0;
		    ct1 = 0.;
		}
		
		ampl[i1]=
		    scan[it1  ][it0  ][i1]*(1.-ct0)*(1.-ct1) +
		    scan[it1+1][it0  ][i1]*(1.-ct0)*ct1 +
		    scan[it1  ][it0+1][i1]*ct0*(1.-ct1) +
		    scan[it1+1][it0+1][i1]*ct0*ct1;
	    }
	} else {
	    for (i1=0; i1 < n1; i1++) {
		pick2[0][i1] = o2+pick2[0][i1]*d2;
		pick2[1][i1] = o3+pick2[1][i1]*d3;
	    }
	}
	
	if (smooth) {
	    /* normalize amplitudes */
	    asum = 0.;
	    for (i1 = 0; i1 < n1; i1++) {
		a = ampl[i1];
		asum += a*a;
	    }
	    asum = sqrtf (asum/n1);
	    for(i1=0; i1 < n1; i1++) {
		ampl[i1] /= asum;
		pick[0][i1] = (o2+pick2[0][i1]*d2-vel2)*ampl[i1];
		pick[1][i1] = (o3+pick2[1][i1]*d3-vel3)*ampl[i1];
	    }
	    
	    sf_divn(pick[0],ampl,pick2[0]);
	    sf_divn(pick[1],ampl,pick2[1]);
	    
	    for(i1=0; i1 < n1; i1++) {
		pick2[0][i1] += vel2;
		pick2[1][i1] += vel3;
	    }
	}
	
	sf_floatwrite(pick2[0],n1*2,pik);
    }

    exit(0);
}
