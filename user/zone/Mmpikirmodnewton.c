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
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <rsf.h>
#include <mpi.h>
#include <omp.h>

#include "mpikirmod.h"
#include "mpikirmod2.h"
#include "mpikirmodnewton.h"
#include "mpikirmodnewton2.h"

#include "mpiricker.h"

int main(int argc, char* argv[]) 
{
    /*For MPI-----------------*/
    MPI_Init(&argc,&argv);
    int size, rank;
    int tstart,tstop;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    /*Timing for the whole algorithm*/
    if(rank == 0) {  
	tstart = MPI_Wtime();    
    }

    /* For newton-------------*/
    int niter, vstatus, order, count, count1, count2, count3;
    double tolerance;
    float **rr, **rd, *updown;
    bool newton, debug, fwdxini;
    velocity2 vn; 
    /*------------------------*/

    int nx, nt, ns, nh, nc, nxc, is, ih, ix, ic;
    float **rfl, **rgd, **crv, **dip, *trace, *trace2;
    float **time, **ampl, **delt, freq, theta, ava, amp, obl;
    float slow, dx, x0, dt, t0, ds, s0, dh, h0, r0;
    const char *type, *type2;
    bool twod, verb, adj, lin, cmp, absoff;
    surface inc, ref;
    velocity vel, vel2;
    ktable ts, tg, **tss, **tgs;
    sf_file data, refl, curv, modl, vti;
	
    sf_init(argc,argv);
	
    if (!sf_getbool("newton",&newton)) newton=true;
    /* To switch between analytical and newton kirmod*/
	
	
    if (!sf_getbool("lin",&lin)) lin=false;
    /* if linear operator */
	
    if (lin) {
	if (!sf_getbool("adj",&adj)) adj=false;
	/* adjoint flag */
		
	if (adj) {
	    modl = sf_input("input");
	    data = sf_output("output");
	} else {
	    data = sf_input("input");
	    modl = sf_output("output");
	}
	curv = sf_input("curv");
    } else {
	adj = false;
	curv = sf_input("input");
	modl = sf_output("output");
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
	
    if (newton) {
	
	if (!sf_getbool("debug",&debug)) debug=false;
	/* debug flag */
		
	if (!sf_getbool("fwdxini",&fwdxini)) fwdxini=false;
	/* use the result of previous iteration to be the xinitial of the next one */
		
	rr = sf_floatalloc2(nx,nc+1);
	rd = sf_floatalloc2(nx,nc+1);
		
	for (count2=0; count2<nc+1; count2++) { /* Generate the reflector input for newton*/
	    for (count1=0; count1<nx; count1++) {
		if (count2==0) {
		    rr[count2][count1] = 0; /* For the surface*/
		    rd[count2][count1] = 0;
		}
		else {
		    rr[count2][count1] = crv[count2-1][count1]; 
		    rd[count2][count1] = dip[count2-1][count1]; 
		}
	    }
	}
		
	updown = sf_floatalloc(nc); /* Fix this if need any other multiple*/
		
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
	else {
		vti = sf_input("aniso"); /* anisotropy*/
		sf_floatread(vn.aniso[0],4*(nc),vti);
	}
	
	if (!sf_getint("niter",&niter)) niter=500;
	/* The number of iterations*/
		
if (!sf_getdouble("tol",&tolerance)) tolerance=0.00001;
	/* Assign a default value for tolerance*/
		
	if (!sf_getint("order",&order)) order=3;/* Interpolation order*/

		
    }
    else {
	/*** Initialize velocity ***/
	vel  = (velocity) sf_alloc(1,sizeof(*vel));
	vel2 = (velocity) sf_alloc(1,sizeof(*vel2));
		
	if (!sf_getfloat("vel",&(vel->v0))) sf_error("Need vel=");
	/*( vel velocity )*/
		
	if (!sf_getfloat("gradx",&(vel->gx))) (vel->gx)=0.;
	/*( gradx horizontal velocity gradient )*/
		
	if (!sf_getfloat("gradz",&(vel->gz))) (vel->gz)=0.;
	/*( gradz vertical velocity gradient )*/
		
	if (!sf_getfloat("velz",&(vel->vz))) (vel->vz)=vel->v0;
	/*( velz vertical velocity for VTI anisotropy )*/
	if (!sf_getfloat("eta",&(vel->n))) (vel->n)=0.;
	/*( eta parameter for VTI anisotropy )*/
		
	type = sf_getstring("type");
	/* type of velocity, 'c': constant, 's': linear sloth, 'v': linear velocity, 'a': VTI anisotropy */
	if (NULL==type) {
	    type= ((vel->gx)==0. && (vel->gz)==0.)?"const":"veloc";
	} else if ((vel->gx)==0. && (vel->gz)==0. && (vel->n)==0.) {
	    type = "const"; 
	} else if ('s'==type[0]) {
	    /* linear slowness squared */
			
	    slow = 1./((vel->v0)*(vel->v0));
	    /* slowness squared */
	    (vel->gx) *= -2.*slow/(vel->v0);
	    (vel->gz) *= -2.*slow/(vel->v0);
	    (vel->v0) = slow;     
	} else if ('v' != type[0] && 'a' != type[0]) {
	    sf_error("Unknown type=%s",type);
	}
		
	if (!sf_getbool("twod",&twod)) twod=false;
	/* 2-D or 2.5-D */
		
	if (!sf_getfloat("refx",&(vel->x0))) (vel->x0)=x0;
	/*( refx reference x-coordinate for velocity )*/
	if (!sf_getfloat("refz",&(vel->z0))) (vel->z0)=0.;
	/*( refz reference z-coordinate for velocity )*/
		
	if (!sf_getfloat("vel2",&(vel2->v0))) (vel2->v0)=(vel->v0);
	/*( vel2 converted velocity )*/
		
	if (!sf_getfloat("gradx2",&(vel2->gx))) (vel2->gx)=(vel->gx);
	/*( gradx2 converted velocity, horizontal gradient )*/
	if (!sf_getfloat("gradz2",&(vel2->gz))) (vel2->gz)=(vel->gz);
	/*( gradz2 converted velocity, vertical gradient )*/
		
		
	type2 = sf_getstring("type2");
	/* type of velocity for the converted (receiver side) branch */
	if (NULL==type2) {	
	    type2=type;
	} else if ((vel2->gx)==0. && (vel2->gz)==0. && (vel2->n)==0.) {
	    type2 = "const"; 
	} else if ('s'==type2[0]) {
	    /* linear slowness squared */
			
	    slow = 1./((vel2->v0)*(vel2->v0));
	    /* slowness squared */
	    (vel2->gx) *= -slow/(vel2->v0);
	    (vel2->gz) *= -slow/(vel2->v0);
	    (vel2->v0) = slow;     
	} else if ('v' != type2[0] && 'a' != type2[0]) {
	    sf_error("Unknown type=%s",type2);
	}
		
	if (!sf_getfloat("refx2",&(vel2->x0))) (vel2->x0)=(vel->x0);
	if (!sf_getfloat("refz2",&(vel2->z0))) (vel2->z0)=(vel->z0);
	/* reference coordinates for converted velocity */
    }

    /*** Allocate space ***/
	
    if (!sf_getbool("cmp",&cmp)) cmp=false;
    /* compute CMP instead of shot gathers */
    
	
    if (newton) {
	inc = kirmodnewton2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc, cmp, absoff);
	ref = inc;
    }
    else {
	if (cmp && !adj) sf_putint(modl,"CDPtype",1);
		
	inc = kirmod2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc, cmp, absoff);
	if (strcmp(type,type2) ||
	    (vel2->v0) != (vel->v0) || 
	    (vel2->gz) != (vel->gz) ||
	    (vel2->gx) != (vel->gx) ||
	    (vel2->z0) != (vel->z0) ||
	    (vel2->x0) != (vel->x0) ) {
	    ref = kirmod2_init(ns, s0, ds, nh, h0, dh, nx, x0, dx, nc, cmp, absoff);
	} else {
	    ref = inc;
	}
    }
    
    /*** Initialize stretch ***/
    sf_aastretch_init (false, nt, t0, dt, nxc);
	
    time = sf_floatalloc2(nx,nc);
    ampl = sf_floatalloc2(nx,nc);
    delt = sf_floatalloc2(nx,nc);
	
    if (!sf_getfloat("freq",&freq)) freq=0.2/dt;
    /* peak frequency for Ricker wavelet */
    ricker_init(nt*2,freq*dt,2);
	
	
    if (newton) {
	/* Initialize parameters for newton*/
	kirmodnewton_init(rr, rd, updown, x0, dx, nx, nc-1, order, nc+1, vstatus, vn.xref, vn.zref, vn.v, vn.gx, vn.gz, vn.aniso);
	
	/*** Compute traveltime table ***/
	kirmodnewton2_table(inc, debug /* Debug Newton */, fwdxini,  niter, tolerance, size, rank);
    }
    else {
	/*** Compute traveltime table ***/
		
	kirmod2_table (inc, vel, type[0], twod, crv, dip);
	if (ref != inc) kirmod2_table (ref, vel2, type2[0], twod, crv, dip);
    }

if(rank == 0){ /* Execute serially using the master only*/

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

/*#ifdef _OPENMP*/
/*#pragma omp parallel for default(none) collapse(2) \*/
/*	private(ix,ic,ts,tg) shared(nc,nx,newton,inc,ref,ih,is,time,delt,tss,tgs,dx,nh,cmp)*/
/*#endif*/
	for (ic=0; ic < nc; ic++) {
		for (ix=0; ix < nx; ix++) {
		if (cmp) {
						
			if (newton) {
				ts = kirmodnewton2_map(inc,is,2*ih,  ix,ic);
				tg = kirmodnewton2_map(ref,is,2*ih+1,ix,ic);
			}
			else {
				ts = kirmod2_map(inc,is,2*ih,  ix,ic);
				tg = kirmod2_map(ref,is,2*ih+1,ix,ic);
			}
						
		} else {
						
			if (newton) {
				ts = kirmodnewton2_map(inc,is,nh,ix,ic);
				tg = kirmodnewton2_map(ref,is,ih,ix,ic);
			}
			else {
				ts = kirmod2_map(inc,is,nh,ix,ic);
				tg = kirmod2_map(ref,is,ih,ix,ic);
			}
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
	    }

/*#ifdef _OPENMP*/
/*#pragma omp parallel for default(none) collapse(2) \*/
/*	private(ix,ic,ts,tg,obl,amp,theta,ava) shared(nc,nx,lin,adj,rfl,ampl,dx,ref,inc,tss,tgs,rgd)*/
/*#endif	*/
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
	}
			
	    if (!adj) {
		sf_aastretch_lop (false,false,nxc,nt,ampl[0],trace);
				
		/* convolve with Ricker wavelet */
		sf_freqfilt_lop(false,false,nt,nt,trace,trace2);
				
		sf_floatwrite(trace2,nt,modl); 
	    }
	}
    }
    sf_warning(".");
	
    if (lin && adj) sf_floatwrite(rfl[0],nxc,data);

    /*For MPI*/
    tstop = MPI_Wtime();
    sf_warning("Total time %d \n",tstop-tstart);
    MPI_Finalize();
}
else if (rank != 0)  MPI_Finalize();
    
    return 0;
}

