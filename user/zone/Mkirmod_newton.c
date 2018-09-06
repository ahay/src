/* Kirchhoff 2-D/2.5-D modeling in layered media with bending ray tracing.  */
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

#include <float.h>
#include <math.h>
#include <rsf.h>
#include <time.h>

#include "kirmodnewton.h"
#include "kirmodnewton2.h"

#include "ricker.h"

int main(int argc, char* argv[]) 
{
    /*Timing*/
    clock_t tstart, tstop;
    double timespend;

    /* For newton-------------*/
    int niter, vstatus, order, count, count1, count2, count3;
    double tolerance;
    float **rr=NULL, **rd=NULL;
    int *updown=NULL;
    bool debug, fwdxini;
    velocity2 vn; 
    /*------------------------*/

    int nx, nt, ns, nh, nc, nxc, is, ih, ix, ic, minix;
    float **rfl, **rgd, **crv, **dip, *trace, *trace2;
    float **time, **ampl, **delt, freq, theta, ava, amp, obl;
    float dx, x0, dt, t0, ds, s0, dh, h0, r0, mint;
    bool verb, adj, lin, cmp, absoff;
    surface inc, ref;
    ktable ts, tg, **tss, **tgs;
    sf_file data, refl, curv, modl, vti, picks = NULL, slopes = NULL;

    tstart = clock();
	
    sf_init(argc,argv);	
	
    if (!sf_getbool("lin",&lin)) lin=false;
    /* if linear operator */
	
    if (lin) {
	if (!sf_getbool("adj",&adj)) adj=false;
	/* adjoint flag */
		
	if (adj) {
	    modl = sf_input("in");
	    data = sf_output("out");
	} else {
	    data = sf_input("in");
	    modl = sf_output("out");
	}
	curv = sf_input("curv");
    } else {
	adj = false;
	curv = sf_input("in");
	modl = sf_output("out");
	data = NULL;
    }
    
    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint  (curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histint  (curv,"n2",&nc)) nc=1; /* number of reflectors */
    nxc = nx*nc;
    
    if (!sf_getbool("absoff",&absoff)) absoff=false;
    /* y - h0 is not in shot coordinate system */
	
    /*** Initialize trace ***/
    
    if (adj) {
	if (!sf_histint(modl,"n1",&nt)) sf_error("Need nt=");
	if (!sf_histfloat(modl,"d1",&dt)) dt=0.004;
	if (!sf_histfloat(modl,"o1",&t0)) t0=0.;
		
		
	/*** Initialize shots ***/
		
	if (!sf_histint(modl,"n3",&ns)) ns=nx;
	if (!sf_histfloat(modl,"o3",&s0)) s0=x0;
	if (!sf_histfloat(modl,"d3",&ds)) ds=dx;
		
	/*** Initialize offsets ***/
		
	if (!sf_histint  (modl,"n2",&nh)) nh=nx;
	if (!sf_histfloat(modl,"o2",&h0)) h0=0.;
	if (!sf_histfloat(modl,"d2",&dh)) dh=dx;
		
		
	sf_putint(data,"n1",nx);
	sf_putfloat(data,"d1",dx);
	sf_putfloat(data,"o1",x0);
	sf_putint(data,"n2",nc);
	sf_putint(data,"n3",1);
    } else {
	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples */
	if (!sf_getfloat("dt",&dt)) dt=0.004;
	/* time sampling */
	if (!sf_getfloat("t0",&t0)) t0=0.;
	/* time origin */
		
	sf_putint  (modl,"n1",nt);
	sf_putfloat(modl,"d1",dt);
	sf_putfloat(modl,"o1",t0);
	sf_putstring(modl,"label1","Time");
	sf_putstring(modl,"unit1","s");
		
	/*** Initialize shots ***/
		
	if (!sf_getint("ns",&ns)) ns=nx;
	/* number of shots (midpoints if cmp=y) */
	if (!sf_getfloat("s0",&s0)) s0=x0;
	/* first shot (midpoint if cmp=y) */
	if (!sf_getfloat("ds",&ds)) ds=dx;
	/* shot/midpoint increment */
		
	sf_putfloat(modl,"o3",s0);
	sf_putfloat(modl,"d3",ds);
	sf_putint(modl,"n3",ns);
		
	/*** Initialize offsets ***/
		
	if (!sf_getint  ("nh",&nh)) nh=nx;
	/* number of offsets */
	if (!sf_getfloat("h0",&h0)) h0=0.;
	/* first offset */
	if (!sf_getfloat("dh",&dh)) dh=dx;
	/* offset increment */
		
	sf_putint  (modl,"n2",nh);
	sf_putfloat(modl,"o2",h0);
	sf_putfloat(modl,"d2",dh);
    }
	
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
	
    trace = sf_floatalloc(nt);
    trace2 = sf_floatalloc(nt);
	
    /*** Initialize reflector ***/
	
    crv = sf_floatalloc2(nx,nc);
    rfl = sf_floatalloc2(nx,nc);
    rgd = sf_floatalloc2(nx,nc);
    dip = sf_floatalloc2(nx,nc);
    
    sf_floatread(crv[0],nxc,curv);
    
    if (!lin) {
	/* reflectivity (A) */
	if (NULL != sf_getstring("refl")) {
	    refl = sf_input("refl");
	    sf_floatread(rfl[0],nxc,refl);
	    sf_fileclose(refl);
	} else {
	    if (!sf_getfloat("r0",&r0)) r0=1.;
	    /* normal reflectivity (if constant) */
	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    rfl[ic][ix] = r0;
		}
	    }
	}
    }
	
    if (NULL != sf_getstring("rgrad")) {
	/* AVO gradient file (B/A) */
	refl = sf_input("rgrad");
	sf_floatread(rgd[0],nxc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (ix=0; ix < nx; ix++) {
		rgd[ic][ix] = 0.;
	    }
	}
    }
	
    if (NULL != sf_getstring("dip")) {
	/* reflector dip file */
	refl = sf_input("dip");
	sf_floatread(dip[0],nxc,refl);
	sf_fileclose(refl);
    } else {
	for (ic=0; ic < nc; ic++) {
	    for (ix=0; ix < nx; ix++) {
		dip[ic][ix] = 0.;
	    }
	}
    }
	
    if (!lin && !adj) {
	/* reflectivity (A) */
	if (NULL != sf_getstring("refl")) {
	    refl = sf_input("refl");
	    sf_floatread(rfl[0],nxc,refl);
	    sf_fileclose(refl);
	} else {
	    if (!sf_getfloat("r0",&r0)) r0=1.;
	    /* normal reflectivity (if constant) */
	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    rfl[ic][ix] = r0;
		}
	    }
	}
    }
	
    if (!lin && !adj) {
	if (NULL != sf_getstring("picks")) {
	    picks = sf_output("picks");
	    /* Output traveltime picks */
	    sf_putint  (picks,"n1",nc);
	    sf_putfloat(picks,"o1",0);
	    sf_putfloat(picks,"d1",1.0);
	    sf_putint  (picks,"n2",nh);
	    sf_putfloat(picks,"o2",h0);
	    sf_putfloat(picks,"d2",dh);
	    sf_putint  (picks,"n3",ns);
	    sf_putfloat(picks,"o3",s0);
	    sf_putfloat(picks,"d3",ds);
	}
	if (NULL != sf_getstring("slopes")) {
	    slopes = sf_output("slopes");
	    /* Output receiver slopes */
	    sf_putint  (slopes,"n1",nc);
	    sf_putfloat(slopes,"o1",0.0f);
	    sf_putfloat(slopes,"d1",1.0f);
	    sf_putint  (slopes,"n2",nh);
	    sf_putfloat(slopes,"o2",h0);
	    sf_putfloat(slopes,"d2",dh);
	    sf_putint  (slopes,"n3",ns);
	    sf_putfloat(slopes,"o3",s0);
	    sf_putfloat(slopes,"d3",ds);
	}
    }
	
    if (!sf_getbool("debug",&debug)) debug=false;
    /* debug flag */
		
    if (!sf_getbool("fwdxini",&fwdxini)) fwdxini=false;
    /* use the result of previous iteration to be the xinitial of the next one */
		
    rr = sf_floatalloc2(nx,nc+1);
    rd = sf_floatalloc2(nx,nc+1);
		
    for (count2=0; count2<nc+1; count2++) { /* Generate the reflector input for newton*/
	for (count1=0; count1<nx; count1++) {
	    if (count2==0) {
		rr[count2][count1] = 0; /* For generating the surface*/
		rd[count2][count1] = 0;
	    }
	    else {
		rr[count2][count1] = crv[count2-1][count1]; 
		rd[count2][count1] = dip[count2-1][count1]; 
	    }
	}
    }
		
    updown = sf_intalloc(nc); /* Fix this if need any other multiple*/
		
    for (count3=0; count3<nc; count3++) {
	updown[count3] = count3+1;
    }
		
    vn.v = sf_floatalloc(nc);
    vn.xref = sf_floatalloc(nc);
    vn.zref = sf_floatalloc(nc);
    vn.gx = sf_floatalloc(nc);
    vn.gz = sf_floatalloc(nc);
    vn.aniso = sf_floatalloc2(4,nc);
	
	
    if (!sf_getint("vstatus",&vstatus)) sf_error("Please enter the status of velocity (0 for constant v,1 for gradient v, and 2 for VTI)");
    /* Velocity status (0 for constant v,1 for gradient v, and 2 for vti)*/
	
    if (vstatus != 2) {
	if (!sf_getfloats("velocity",vn.v,nc)) sf_error("Please enter the velocity array [nc]");
	/* Assign velocity km/s*/
	if (vstatus == 1) {
	    if (!sf_getfloats("xgradient",vn.gx,nc)) {
		for (count=0; count<nc; count++) {
		    vn.gx[count] = 0;
		}
	    }
	    /* Assign x-gradient*/
				
	    if (!sf_getfloats("zgradient",vn.gz,nc)) { 
		for (count=0; count<nc; count++) {
		    vn.gz[count] = 0;
		}
	    }
	    /* Assign z-gradient */
				
	    if (!sf_getfloats("xref",vn.xref,nc))  sf_error("Please enter the x-reference points array [nc]");
	    /* Assign x-reference point*/
				
	    if (!sf_getfloats("zref",vn.zref,nc)) sf_error("Please enter the z-reference points array [nc]");
	    /* Assign z-reference point*/
	}
    }
    else {
	vti = sf_input("aniso"); /* anisotropy*/
	sf_floatread(vn.aniso[0],4*(nc),vti);
    }
	
    if (!sf_getint("niter",&niter)) niter=500;
    /* The number of iterations*/
		
    if (!sf_getdouble("tol",&tolerance)) tolerance=0.00001;
    /* Assign a default value for tolerance*/
		
    if (!sf_getint("order",&order)) order=3;/* Interpolation order*/

		
    /*** Allocate space ***/
	
    if (!sf_getbool("cmp",&cmp)) cmp=false;
    /* compute CMP instead of shot gathers */
    
    inc = kirmodnewton2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc, cmp, absoff);
    ref = inc;
 
    /*** Initialize stretch ***/
    sf_aastretch_init (false, nt, t0, dt, nxc);
	
    time = sf_floatalloc2(nx,nc);
    ampl = sf_floatalloc2(nx,nc);
    delt = sf_floatalloc2(nx,nc);
	
    if (!sf_getfloat("freq",&freq)) freq=0.2/dt;
    /* peak frequency for Ricker wavelet */
    ricker_init(nt*2,freq*dt,2);

    /* Initialize parameters for newton*/
    kirmodnewton_init(rr, rd, updown, x0, dx, nx, nc-1, order, nc+1, vstatus, vn.xref, vn.zref, vn.v, vn.gx, vn.gz, vn.aniso);
    /*** Compute traveltime table ***/
    kirmodnewton2_table(inc, debug /* Debug Newton */, fwdxini,  niter, tolerance);

    if (lin) {
	if (adj) {
	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    rfl[ic][ix] = 0.0;
		}
	    }
	} else {
	    sf_floatread(rfl[0],nxc,data);
	}
    }
	
    tss = (ktable**) sf_alloc(nc,sizeof(ktable*));
    tgs = (ktable**) sf_alloc(nc,sizeof(ktable*));
    for (ic=0; ic < nc; ic++) {
	tss[ic] = (ktable*) sf_alloc(nx,sizeof(ktable));
	tgs[ic] = (ktable*) sf_alloc(nx,sizeof(ktable));
    }

    /*** Main loop ***/
    for (is=0; is < ns; is++) {
	if (verb) sf_warning("%s %d of %d;",cmp?"cmp":"shot",is+1,ns);
		
	for (ih=0; ih < nh; ih++) {
			
	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    if (cmp) {
			ts = kirmodnewton2_map(inc,is,2*ih,  ix,ic);
			tg = kirmodnewton2_map(ref,is,2*ih+1,ix,ic);
		    } else {
			ts = kirmodnewton2_map(inc,is,nh,ix,ic);
			tg = kirmodnewton2_map(ref,is,ih,ix,ic);
		    }
					
		    time[ic][ix] = ts->t + tg->t;
		    delt[ic][ix] = fabsf(ts->tx+tg->tx)*dx;
					
		    tss[ic][ix] = ts;
		    tgs[ic][ix] = tg;
		}
	    }
			
	    sf_aastretch_define (time[0],delt[0],NULL);
			
	    if (adj) {
		sf_floatread(trace2,nt,modl);
				
		/* convolve with Ricker wavelet */
		sf_freqfilt_lop(true,false,nt,nt,trace,trace2);
				
		sf_aastretch_lop (true,false,nxc,nt,ampl[0],trace); 
				
                /* aastretch_lop (true,false,nxc,nt,ampl[0],trace2); */
	    }
			
	    for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		    ts = tss[ic][ix];
		    tg = tgs[ic][ix];
					
		    obl = 0.5*(ts->tn + tg->tn);
		    amp = ts->a * tg->a * sqrtf(ts->ar + tg->ar) + FLT_EPSILON;
					
		    if (lin) {
			if (adj) {
			    rfl[ic][ix] += ampl[ic][ix]*obl*dx/amp;
			} else {
			    ampl[ic][ix] = rfl[ic][ix]*obl*dx/amp;
			}
		    } else {
			theta = 0.5*(SF_SIG(tg->tx)*tg->an - 
				     SF_SIG(ts->tx)*ts->an);
			theta = sinf(theta);
						
			ava = rfl[ic][ix]+rgd[ic][ix]*theta*theta;
			if (ref != inc) ava *= theta;
						
			ampl[ic][ix] = ava*obl*dx/amp;
		    }
		}
		/* Pick traveltime and/or receiver slope */
		if (picks || slopes) {
		    mint = SF_HUGE;
		    minix = 0;
		    for (ix=0; ix < nx; ix++) {
			if (time[ic][ix] < mint) {
			    mint = time[ic][ix];
			    minix = ix;
			}
		    }
		    if (picks) {
			amp = tss[ic][minix]->t + tgs[ic][minix]->t;
			sf_floatwrite(&amp,1,picks);
		    }
		    if (slopes) {
			amp = fabsf(tgs[ic][minix]->tx);
			sf_floatwrite(&amp,1,slopes);
		    }
		}
	    }
			
	    if (!adj) {
		sf_aastretch_lop (false,false,nxc,nt,ampl[0],trace);
				
		/* convolve with Ricker wavelet */
		sf_freqfilt_lop(false,false,nt,nt,trace,trace2);
				
		sf_floatwrite(trace2,nt,modl); 
				
                /* sf_floatwrite(trace,nt,modl); */
	    }
	}
    }
    sf_warning(".");
	
    if (lin && adj) sf_floatwrite(rfl[0],nxc,data);
	
    tstop = clock();
    timespend = (double)(tstop-tstart)/CLOCKS_PER_SEC;
	
    sf_warning("Total computational time %g",timespend);
	
    exit(0);
}

