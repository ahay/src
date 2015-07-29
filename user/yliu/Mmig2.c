/* 2-D Prestack Kirchhoff time migration with antialiasing. 
The axes in the input are {time,midpoint,offset}
The axes in the offset are {1,midpoint,offset}
The axes in the output are {time,midpoint}
The axes in the "image gather" are {time,midpoint,offset}
*/
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

static float pick(float ti, float deltat, 
		  const float *trace,
		  int nt, float dt, float t0)
/* pick a traveltime sample from a trace */
{
    int it, itm, itp;
    float ft, tm, tp, ftm, ftp, imp;

    ft = (ti-t0)/dt; it = floorf(ft); ft -= it; 
    if ( it < 0 || it >= nt-1) return 0.0;
 
    tm = ti-deltat-dt;
    ftm = (tm-t0)/dt; itm = floorf(ftm); ftm -= itm; 
    if (itm < 0) return 0.0;
                 
    tp = ti+deltat+dt;
    ftp = (tp-t0)/dt; itp = floorf(ftp); ftp -= itp; 
    if (itp >= nt-1) return 0.0;

    imp = dt/(dt+tp-tm);
    imp *= imp;
		
    return imp*(
	2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
	(1.-ftm)*trace[itm] - ftm*trace[itm+1]    - 
	(1.-ftp)*trace[itp] - ftp*trace[itp+1]);    
}

int main(int argc, char* argv[])
{
    int nt, nx, nh, nh2, ix, ih, iy, i, nn, it, **fold, apt;
    float *trace, **out, **v, rho, **outd, *pp, *off;
    float h, x, t, h0, dh, dx, ti, tx, t0, t1, t2, dt, vi, aal, angle;
    sf_file inp, mig, vel, gather, offset;
    bool half, verb;

    sf_init (argc,argv);
    inp = sf_input("in");
    vel = sf_input("vel");
    mig = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1=");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2=");
     if (!sf_histint(inp,"n3",&nh)) sf_error("No n3=");

    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1=");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1=");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2=");

    sf_putint(mig,"n3",1);
    sf_putfloat(mig,"d3",1.);
    sf_putfloat(mig,"o3",0.);
    
    if (NULL != sf_getstring("gather")) {
	gather = sf_output("gather");
    } else {
	gather = NULL;
    }

    if (!sf_getfloat("antialias",&aal)) aal = 1.0;
    /* antialiasing */

    if (!sf_getint("apt",&apt)) apt = nx;
    /* integral aperture */

    if (!sf_getfloat("angle",&angle)) angle = 90.0;
    /* angle aperture */

    angle = fabsf(tanf(angle*SF_PI/180.0));

    if (!sf_getbool("half",&half)) half = true;
    /* if y, the third axis is half-offset instead of full offset */

    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */

    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	nh2 = sf_filesize(offset);

	if (nh2 != nh*nx) sf_error("Wrong dimensions in offset");

	off = sf_floatalloc(nh2);	
	sf_floatread (off,nh2,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(inp,"o3",&h0)) sf_error("No o3=");
	if (!sf_histfloat(inp,"d3",&dh)) sf_error("No d3=");
	
	if (!half) dh *= 0.5;

	off = sf_floatalloc(nh*nx);
	for (ix = 0; ix < nx; ix++) {
	    for (ih = 0; ih < nh; ih++) {
		off[ih*nx+ix] = h0 + ih*dh; 
	    }
	}
	offset = NULL;
    }

    v = sf_floatalloc2(nt,nx);
    sf_floatread(v[0],nt*nx,vel);

    trace = sf_floatalloc(nt);
    out = sf_floatalloc2(nt,nx);
    outd = sf_floatalloc2(nt,nx);
    fold = sf_intalloc2(nt,nx);

    nn = 2*kiss_fft_next_fast_size((nt+1)/2);
    pp = sf_floatalloc(nn);

    sf_halfint_init (true, nn, rho);

    for (i=0; i < nt*nx; i++) {
	outd[0][i] = 0.;  
	fold[0][i] = 0;  
    }

    for (ih=0; ih < nh; ih++) {
        if (verb) sf_warning("offset %d of %d;",ih+1,nh);
	for (i=0; i < nt*nx; i++) {
	    out[0][i] = 0.;
	}

	for (iy=0; iy < nx; iy++) { 
	    sf_floatread (trace,nt,inp);
	    sf_doubint(true, nt,trace);
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

		    out[ix][it] += 
			pick(ti,fabsf(tx*dx*aal),trace,nt,dt,t0);
		} 
	    } 
	} 

	for (iy=0; iy < nx; iy++) {
	    for (it=0; it < nt; it++) {
		pp[it] = out[iy][it];
	    }
	    for (it=nt; it < nn; it++) {
		pp[it] = 0.;
	    }
	    sf_halfint (true, pp);
	    for (it=0; it < nt; it++) {
		out[iy][it] = pp[it];
	    }
	}

	if (NULL != gather) sf_floatwrite(out[0],nt*nx,gather);

	for (iy=0; iy < nx; iy++) {
	    for (it=0; it < nt; it++) {
		outd[iy][it] += out[iy][it];
		if (0.!=out[iy][it]) fold[iy][it]++;
	    }
	}
    }
    if (verb) sf_warning(".");

    for (i=0; i < nt*nx; i++) {
	outd[0][i] = outd[0][i]/(fold[0][i]+FLT_EPSILON);  
    }
    sf_floatwrite(outd[0],nt*nx,mig); 

    exit(0);
}

/* 	$Id: Mmig2.c 13985 2015-03-26 13:56:59Z sfomel $	 */
