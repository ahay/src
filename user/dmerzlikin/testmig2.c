/* Kirchhoff zero-offset modeling/migration antialiased by parameterization */
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
//////*^*/

#include "testmig2.h"

void pick(bool adj, float ti, float deltat, 
		 float *trace, float *out, int i,
		 int nt, float dt, float t0)
/*< pick a traveltime sample from a trace >*/
{
    int it, itm, itp;
    float ft, tm, tp, ftm, ftp, imp;

    ft = (ti-t0)/dt; it = floorf(ft); ft -= it; 
    if ( it < 0 || it >= nt-1) return;
 
    tm = ti-deltat-dt;
    ftm = (tm-t0)/dt; itm = floorf(ftm); ftm -= itm; 
    if (itm < 0) return;
                 
    tp = ti+deltat+dt;
    ftp = (tp-t0)/dt; itp = floorf(ftp); ftp -= itp; 
    if (itp >= nt-1) return;

    imp = dt/(dt+tp-tm);
    imp *= imp;

    if (adj) {
	out[i] += imp*(
	    2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
	    (1.-ftm)*trace[itm] - ftm*trace[itm+1]    - 
	    (1.-ftp)*trace[itp] - ftp*trace[itp+1]);    
    } else {
	trace[it] += 2.*(1-ft)*imp*out[i];
	trace[it+1] += 2.*ft*imp*out[i];
	trace[itm] -= (1.-ftm)*imp*out[i];
	trace[itm+1] -= ftm*imp*out[i];
	trace[itp] -= (1.-ftp)*imp*out[i];
	trace[itp+1] -= ftp*imp*out[i];
    }
}

void mig2_lop (bool adj, bool half, bool verb, bool normalize,
		   int nt, int nx, int nh, int **fold, int apt,
		   float *trace_temp, float **v, float rho, float *stack, float *off,
		   float h0, float dh, float dx, float t0, float dt, float aal, float angle)
/*< Apply >*/
{
    
    int ix, ih, iy, i, it, r, nn; // loop variables
    float *trace, h, x, t, tx, t1, t2, vi, ti;
    float **image, *pp;

    image = sf_floatalloc2(nt,nx);

    trace = sf_floatalloc(nt);

    if (normalize) {
	fold = sf_intalloc2(nt,nx);
    } else {
	fold = NULL;
    }

    nn = 2*kiss_fft_next_fast_size((nt+1)/2);
    pp = sf_floatalloc(nn);

    sf_halfint_init (true, nn, rho);

    /*if (adj) {
	for (i=0; i < nt*nx; i++) {
	    stack[0][i] = 0.;  
	}
    } else {
	sf_floatread(stack[0],nt*nx,inp); 
    }*/

    if (NULL != fold) {
	for (i=0; i < nt*nx; i++) {
	    fold[0][i] = 0;  
	}
    }

    if(nh!=1){
       sf_warning("do not currently work with multiple offsets!!!!!");
    }

    for (ih=0; ih < nh; ih++) {
        if (verb) sf_warning("offset %d of %d;",ih+1,nh);

	if (adj) {
	    for (i=0; i < nt*nx; i++) {
		image[0][i] = 0.;
	    }
	} else {
	    for (iy=0; iy < nx; iy++) {
		for (it=0; it < nt; it++) {
		    image[iy][it] = stack[iy*nt + it];
		}
	    }
	}

	if (!adj) {
	    for (iy=0; iy < nx; iy++) {
		for (it=0; it < nt; it++) {
		    pp[it] = image[iy][it];
		}
		for (it=nt; it < nn; it++) {
		    pp[it] = 0.;
		}
		sf_halfint (false, pp);
		for (it=0; it < nt; it++) {
		    image[iy][it] = pp[it];
		}
	    }
	}

	for (iy=0; iy < nx; iy++) { 
	    if (adj) {
		//sf_floatread (trace,nt,inp);

			for(r=0;r<nt;r++){
				trace[r] = trace_temp[iy*nt + r];
			}

		sf_doubint(true, nt,trace);
	    } else {
		for (it=0; it < nt; it++) {
		    trace[it]=0.0f;
		}
	    }

	    h = fabsf(off[ih*nx+iy]);

	    for (ix=0; ix < nx; ix++) { 
	        x = (ix-iy)*dx;
		if (SF_ABS(ix-iy) > apt) continue;

		for (it=0; it < nt; it++) {
		    t = t0 + it*dt;  
		    vi = v[ix][it];

		    if (fabsf(x) > angle*vi*t) continue;

		    /* hypot(a,b) = sqrt(a*a+b*b) */
		    t1 = hypotf(0.5*t,(x-h)/vi);
		    t2 = hypotf(0.5*t,(x+h)/vi);
		    ti = t1+t2;

		    /* tx = |dt/dx| */
		    tx = fabsf(x-h)/(vi*vi*(t1+dt))+
		         fabsf(x+h)/(vi*vi*(t2+dt));

		    pick(adj,ti,fabsf(tx*dx*aal),trace,image[ix],it,nt,dt,t0);
		} 
	    } 

	    if (!adj) {
		sf_doubint(true, nt,trace);
		//sf_floatwrite (trace,nt,out);

			for(r=0;r<nt;r++){
				trace_temp[iy*nt + r] = trace[r];
			}

	    }
	} 

	if (adj) {
	    for (iy=0; iy < nx; iy++) {
		for (it=0; it < nt; it++) {
		    pp[it] = image[iy][it];
		}
		for (it=nt; it < nn; it++) {
		    pp[it] = 0.;
		}
		sf_halfint (true, pp);
		for (it=0; it < nt; it++) {
		    image[iy][it] = pp[it];
		}
	    }

	    //if (NULL != gather) sf_floatwrite(image[0],nt*nx,gather);

	    for (iy=0; iy < nx; iy++) {
		for (it=0; it < nt; it++) {
		    stack[iy*nt + it] += image[iy][it];
		    //stack[iy][it] = trace_temp[iy*nt + it];
		    if (NULL != fold && 0.!=image[iy][it]) fold[iy][it]++; 
		}
	    }
	}
    }
    if (verb) sf_warning(".");

    if (NULL != fold) {
	for (i=0; i < nt*nx; i++) {
	    stack[i] /= (fold[0][i]+FLT_EPSILON);  
	}
    }
    
    //if (adj) sf_floatwrite(stack[0],nt*nx,out);

}

