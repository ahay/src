/* 2-D seislet transform using frequency-wavenumber offset-azimuth continuation 
Forward transform (adj=n inv=y/n) m=T[d]
Inverse transform (adj=y inv=y)   d=T^(-1)[d]
Adjoint transform (adj=y inv=n)   d=T'[d]
*/
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

#include "fkoclet3.h"

int main(int argc, char *argv[])
{
    int nkx, nky, ih, nh, ia, na, iw, nw, i6, n6, ikx, iky;
    float kx0, ky0, dkx, dky, h0, hi, dh, da, w0, a0, ai, dw, w, kx, ky, eps, maxe; 
    bool inv, verb, adj, dwt, amp;

    char *type, *dir;
    sf_complex *pp, *qq, *pp1, *pp2, *qq1, *qq2;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nh))  sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&na))  sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nkx)) sf_error("No n3= in input");
    if (!sf_histint(in,"n4",&nky)) sf_error("No n4= in input");
    if (!sf_histint(in,"n5",&nw))  sf_error("No n5= in input");

    if (!sf_histfloat(in,"d1",&dh))  sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&da))  sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dkx)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"d4",&dky)) sf_error("No d4= in input");
    if (!sf_histfloat(in,"d5",&dw))  sf_error("No d5= in input");

    if (!sf_histfloat(in,"o1",&h0))  sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&a0))  sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&kx0)) sf_error("No o3= in input");
    if (!sf_histfloat(in,"o4",&ky0)) sf_error("No o4= in input");
    if (!sf_histfloat(in,"o5",&w0))  sf_error("No o5= in input");

    a0 *= SF_PI/180.;
    da *= SF_PI/180.;

    n6 = sf_leftsize(in,5);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, do adjoint transform */

    if (!sf_getbool("dwt",&dwt)) dwt=false;
    /* if y, do wavelet transform */

    if (!sf_getbool("amp",&amp)) amp=false;
    /* if y, true amplitudes continuation */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getfloat("maxe",&maxe)) maxe=10.;
    /* stability constraint */

    if (NULL == (dir=sf_getstring("dir"))) dir="both";
    /* [azimuth,offset,both] direction, the default is both directions  */

    pp = sf_complexalloc(nh*na);   /* data space */
    pp1 = sf_complexalloc(nh);   /* data space (offset direction) */
    pp2 = sf_complexalloc(na);   /* data space (azimuth direction) */

    qq = sf_complexalloc(nh*na);   /* model space */
    qq1 = sf_complexalloc(nh);   /* model space (offset direction) */
    qq2 = sf_complexalloc(na);   /* model space (azimuth direction) */

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */

    /* loop over n6 */
    for (i6=0; i6 < n6; i6++) {
	for (iw=0; iw < nw; iw++) { /* loop over frequency */
	    if (verb) sf_warning("frequency %d of %d",iw+1,nw);
	    w = w0 + iw*dw;
	    for (iky=0; iky < nky; iky++) { /* loop over wavenumber Y */
		ky = ky0 + iky*dky;

		for (ikx=0; ikx < nkx; ikx++) { /* loop over wavenumber X */
		    kx = kx0 + ikx*dkx;

		    if (adj) {
			sf_complexread(qq,nh*na,in);
		    } else {
			sf_complexread(pp,nh*na,in);
		    } 
		    
		    if (adj) {
			switch (dir[0]) {
			    case 'o':
				/* Inverse Transform (offset direction) */
				for (ia=0; ia < na; ia++) {
				    for (ih=0; ih < nh; ih++) {
					qq1[ih] = qq[ia*nh+ih];
				    }
				    ai = a0+ia*da;
				    fkoclet_init(nh,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(false,false,nh,nh,qq1,pp1,w,kx,ky,h0,dh,ai,0.);
				    for (ih=0; ih < nh; ih++) {
					pp[ia*nh+ih] = pp1[ih];
				    }
				}
				break;
			    case 'a':
				/* Inverse Transform (azimuth direction) */
				for (ih=0; ih < nh; ih++) {
				    for (ia=0; ia < na; ia++) {
					qq2[ia] = qq[ia*nh+ih];
				    }
				    hi = h0+ih*dh;
				    fkoclet_init(na,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(false,false,na,na,qq2,pp2,w,kx,ky,hi,0,a0,da);
				    for (ia=0; ia < na; ia++) {
					pp[ia*nh+ih] = pp2[ia];
				    }
				}
				break;
			    case 'b':
				/* Inverse Transform (offset direction) */
				for (ia=0; ia < na; ia++) {
				    for (ih=0; ih < nh; ih++) {
					qq1[ih] = qq[ia*nh+ih];
				    }
				    ai = a0+ia*da;
				    fkoclet_init(nh,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(false,false,nh,nh,qq1,pp1,w,kx,ky,h0,dh,ai,0.);
				    for (ih=0; ih < nh; ih++) {
					pp[ia*nh+ih] = pp1[ih];
				    }
				}
				/* Inverse Transform (azimuth direction) */
				for (ih=0; ih < nh; ih++) {
				    for (ia=0; ia < na; ia++) {
					qq2[ia] = pp[ia*nh+ih];
				    }
				    hi = h0+ih*dh;
				    fkoclet_init(na,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(false,false,na,na,qq2,pp2,w,kx,ky,hi,0,a0,da);
				    for (ia=0; ia < na; ia++) {
					pp[ia*nh+ih] = pp2[ia];
				    }
				}
				break;
			    default:
				sf_error("Unknown direction \"%s\"",dir);
			}
			sf_complexwrite(pp,nh*na,out);
		    } else {
			switch (dir[0]) {
			    case 'a':
				/* Forward Transform (azimuth direction) */
				for (ih=0; ih < nh; ih++) {
				    for (ia=0; ia < na; ia++) {
					pp2[ia] = pp[ia*nh+ih];
				    }
				    hi = h0+ih*dh;
				    fkoclet_init(na,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(true,false,na,na,qq2,pp2,w,kx,ky,hi,0,a0,da);
				    for (ia=0; ia < na; ia++) {
					qq[ia*nh+ih] = qq2[ia];
				    }
				}
				break;
			    case 'o':
				/* Forward Transform (offset direction) */
				for (ia=0; ia < na; ia++) {
				    for (ih=0; ih < nh; ih++) {
					pp1[ih] = pp[ia*nh+ih];
				    }
				    ai = a0+ia*da;
				    fkoclet_init(nh,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(true,false,nh,nh,qq1,pp1,w,kx,ky,h0,dh,ai,0.);
				    for (ih=0; ih < nh; ih++) {
					qq[ia*nh+ih] = qq1[ih];
				    }
				}
				break;
			    case 'b':
				/* Forward Transform (azimuth direction) */
				for (ih=0; ih < nh; ih++) {
				    for (ia=0; ia < na; ia++) {
					pp2[ia] = pp[ia*nh+ih];
				    }
				    hi = h0+ih*dh;
				    fkoclet_init(na,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(true,false,na,na,qq2,pp2,w,kx,ky,hi,0,a0,da);
				    for (ia=0; ia < na; ia++) {
					qq[ia*nh+ih] = qq2[ia];
				    }
				}
				/* Forward Transform (offset direction) */
				for (ia=0; ia < na; ia++) {
				    for (ih=0; ih < nh; ih++) {
					pp1[ih] = qq[ia*nh+ih];
				    }
				    ai = a0+ia*da;
				    fkoclet_init(nh,inv,false,dwt,amp,eps*eps,type[0],maxe);
				    fkoclet_lop(true,false,nh,nh,qq1,pp1,w,kx,ky,h0,dh,ai,0.);
				    for (ih=0; ih < nh; ih++) {
					qq[ia*nh+ih] = qq1[ih];
				    }
				}
				break;
			    default:
				sf_error("Unknown direction \"%s\"",dir);
			}
			sf_complexwrite(qq,nh*na,out);
		    } 
		}
	    }
	}
    }

    exit(0);
}
/* 	$Id$	 */
