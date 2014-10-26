/* FK-domain Gardner's DMO for regularly sampled 2-D data 

   The input/ouput is (offset,logstretch frequency, midpoint/shot wavenumber).
*/
/*
  Copyright (C) 2013 University of Texas at Austin

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

#include "stretch4.h"

int main(int argc, char* argv[])
{
    bool inv, shot;
    map4 *map;
    int iw, nw, ih, nh, ik, nk, ib, nb, *fold;
    float dw, dh, dk, w0, h0, k0, w, h, k, eps, db, b, sinb, cosb, xi;
    sf_complex *slice, *slice2, *sum, sample, tshift, xshift; 
    float *hstr, *dt, *dx, *ds;
    sf_file in, out;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("inv",&inv)) inv=true;
    /* inversion flag */

    if (!sf_getbool("shot",&shot)) shot=false;
    /* if shot gathers instead of midpoint gathers */

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nk)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dk)) sf_error("No d3= in input");

    if (!sf_histfloat(in,"o1",&h0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&w0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&k0)) sf_error("No o3= in input");

    dk *= 2*SF_PI;
    k0 *= 2*SF_PI;

    dw *= 2*SF_PI;
    w0 *= 2*SF_PI;

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    slice  = sf_complexalloc(nh);
    hstr   = sf_floatalloc(nh);
    slice2 = sf_complexalloc(nh);

    sum  = sf_complexalloc(nh);
    fold = sf_intalloc(nh);

    if (!sf_getint("nb",&nb)) nb=85;   /* number of angles */
    if (!sf_getfloat("db",&db)) db=1;   /* angle increment */
    db *= SF_PI/180;

    if (!sf_getfloat("xi",&xi)) xi=1;   /* continuation paremeter */

    map = (map4*) sf_alloc(nb,sizeof(map4));
    dt = sf_floatalloc(nb);
    dx = sf_floatalloc(nb);
    ds = shot? sf_floatalloc(nb):NULL;

    /* precompute things */
    for (ib=0; ib < nb; ib++) {
	b = (ib+1)*db;
	sinb = sinf(b);
	cosb = sqrtf(1.0f-xi*sinb*sinb);

	dt[ib] = logf(cosb);

	for (ih=0; ih < nh; ih++) {
	    h = h0 + ih*dh;

	    if (inv) {
		hstr[ih] = h*cosb;
	    } else {
		hstr[ih] = h/cosb;
	    }
	}

	if (inv) {
	    dx[ib]=xi*sinb;
	} else {
	    dx[ib]=xi*sinb/cosb;
	}

	if (shot) {
	    if (inv) {
		ds[ib] = 1.0f-cosb;
	    } else {
		ds[ib] = 1.0f-1.0f/cosb;
	    }
	}

	map[ib] = stretch4_init (nh, h0, dh, nh, eps);
	stretch4_define (map[ib],hstr);
    }	
	
    for (ik=0; ik < nk; ik++) {
	k = k0+ik*dk;   
	sf_warning("wavenumber %d of %d;",ik+1,nk);

	for (iw=0; iw < nw; iw++) {
	    w = w0+iw*dw;

	    sf_complexread(slice,nh,in);

	    for (ih=0; ih < nh; ih++) {
		sum[ih] = slice[ih];
		fold[ih] = 1;
	    }
	
	    for (ib=0; ib < nb; ib++) {
		b = w*dt[ib];
		sinb = sinf(b);
		cosb = sqrtf(1.0f-xi*sinb*sinb);

		tshift = sf_cmplx(cosb,-sinb);

		if (inv) {
		    cstretch4_apply  (map[ib],slice,slice2);
		} else {
		    cstretch4_invert (map[ib],slice2,slice);
		}

		for (ih=0; ih < nh; ih++) {
		    h = h0 + ih*dh;

		    if (shot) {
			b = k*h*ds[ib];
			sinb = sinf(b);
			cosb = sqrtf(1.0f-xi*sinb*sinb);

			xshift = tshift*sf_cmplx(cosb,sinb);
		    } else {
			xshift = tshift;
		    }
		    
		    sample = slice2[ih];

		    if (crealf(sample) != 0.0f || 
			cimagf(sample) != 0.0f) {
			sum[ih]  += sample*xshift*2.0f*cosf(k*h*dx[ib]);
			fold[ih] += 2;
		    }
		} /* ih */
	    } /* ib */
	
	    for (ih=0; ih < nh; ih++) {
		if (fold[ih] > 1) sum[ih] /= fold[ih];
	    }
	    
	    sf_complexwrite(sum,nh,out);
	} /* iw */
    } /* ik */
    sf_warning(".");
    
    exit(0);
}
