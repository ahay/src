/* Rice HPCSS reverse time migration. */

/*************************************************************************

Copyright Rice University, 2009.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

/* Modified for distribution with Madagascar */

#include <rsf.h>

#include "wavefun.h"
#include "step.h"

#define CFL 0.7

int main(int argc, char ** argv) {
    
    /* BEGIN DECLARATIONS */
    
    /* input data */
    
    WINFO ri;              /* struct for command line input */	
    
    /* workspace */
    
    float * v;    /* velocity field */
    float * p1;   /* pressure field, current time step */
    float * p0;   /* pressure field, last time step */
    float * rp1;  /* receiver pressure field, current time step */
    float * rp0;  /* receiver pressure field, last time step */
    float * tr;            /* traces */
    float * tmp;           /* swap */
    float * imag;          /* image */ 
    
    int ix, iz, it;        /* counters */
    int isrc;              /* source counter */
    int isx;               /* source location, in units of dx */
    int nxz;               /* number of spatial grid points */
    int nz;                /* local storage for number of depth points */
    int ioff;              /* array offset */
    int ntr;               /* number of traces */
    int nsam;              /* number of samples in shot record */
    float rz, rx, s;       /* precomputed coefficients */
    float vmax, vmin;      /* max and min velocity */ 
    /* float two;             two */
    
    /* END DECLARATIONS */
    
    /* PARSE ARGUMENTS */
    
    sf_init(argc,argv);
    getinputs(false,&ri);
    
    /* ALLOCATE WORKSPACE */
    
    /* compute number of spatial grid points */
    nxz=ri.nx * ri.nz;
    
    /* compute number of traces, samples in each record */
    ntr=ri.igxend-ri.igxbeg+1;
    nsam=ntr*ri.nt;
    
    /* allocate, initialize p0, p1, v, traces, rp0, rp1 */
    p0=sf_floatalloc(nxz);
    p1=sf_floatalloc(nxz);
    rp0=sf_floatalloc(nxz);
    rp1=sf_floatalloc(nxz);
    imag=sf_floatalloc(nxz);
    v=sf_floatalloc(nxz);
    tr=sf_floatalloc(nsam);

    /* INITIALIZE FILE I/O */
    sf_floatread(v,nxz,ri.vfile);

    /* CFL, sanity checks */
    vmax=fgetmax(v,nxz);
    vmin=fgetmin(v,nxz);
    if (vmax*ri.dt>CFL*fmaxf(ri.dx,ri.dz)) {
	sf_warning("ERROR: CFL criterion violated");
	sf_warning("vmax=%e dx=%e dz=%e dt=%e",vmax,ri.dx,ri.dz,ri.dt);
	sf_error("max permitted dt=%e",CFL*fmaxf(ri.dx,ri.dz)/vmax);
    }
    if (vmin<=0.0) sf_error("min velocity nonpositive");

    /* only square of velocity array needed from here on */
    fsquare(v,nxz);

    /* precalculate some coefficients */
    rz=ri.dt*ri.dt/(ri.dz*ri.dz);
    rx=ri.dt*ri.dt/(ri.dx*ri.dx);
    s =2.0*(rz+rx);
    /* two=2.0; */
    nz=ri.nz;

    /* initialize image field */
    fzeros(imag,nxz);
    
    /* SOURCE LOOP */

    isrc=0;
    isx=ri.isxbeg;

    while (isx <= ri.isxend) { 

	/* FORWARD TIME LOOP */

	/* seek to begin of source wavefield file */
	if (NULL != ri.mfile) sf_seek(ri.mfile,0L,SEEK_SET);

	/* compute the source wavefile for one shot */
	fzeros(p0,nxz);
	fzeros(p1,nxz);
 
	/* precondition: before step it, p0 = time (it-1)*dt,
	   p1 = time it*dt.
	   postcondition: after step it, p0 = time it*dt, 
	   p1 = time (it+1)*dt
	*/

	for(it=0;it<ri.nt;it++) {

	    /* construct next time step, overwrite on p0 */
	    step_forward(p0,p1,v,ri.nz,ri.nx,rz,rx,s);
 
	    p0[ri.isz+isx*nz]+=fgetrick(it*ri.dt,ri.freq);

	    /* swap pointers */
	    tmp=p0;
	    p0=p1;
	    p1=tmp;

	    /* optionally write source wavefield at time it to file -
	       note that after pointer swap this is now p0! 
	    */
	    if (NULL != ri.mfile) sf_floatwrite(p0,nxz,ri.mfile);
	}

	/* REVERSE TIME LOOP */
        
	/* receiver wavefield initialization */
	fzeros(rp0,nxz);
	fzeros(rp1,nxz);
	/* imf = 0; */
  
	/* compute receiver wavefield backward in time and image by
	   crosscorrelation at zero time, space lag */
	
	/* read traces for this shot */
	sf_seek(ri.tfile,nsam*isrc*sizeof(float),SEEK_SET);
	sf_floatread(tr,nsam,ri.tfile);

	/* if source wavefields are to be time-stepped backwards,
	   swap pointers so that p0 = time nt*dt 
	*/
	if (NULL != ri.mfile) {
	    tmp=p0;
	    p0=p1;
	    p1=tmp;
	}
       
	for (it=ri.nt-1;it>-1;it--) {
		
	    /* construct next time step (backward) of receiver wavefield,
	       overwrite on rp0.  do same to source wavefield if not saved
	       to disk
	       
	       precondition: before step it, rp0 = time (it+1)*dt,
	       rp1 = time it*dt.
	       postcondition: after step it, rp0 = time it*dt, 
	       rp1 = time (it-1)*dt
	       same for p0, p1
	    */       
	    if (NULL != ri.mfile) {
		step_forward(rp0,rp1,v,ri.nz,ri.nx,rz,rx,s);

		/* swap pointers */
		tmp=rp0;
		rp0=rp1;
		rp1=tmp;
		
	    } else {
		/* going backwards in time, SUBTRACT source BEFORE step */
		p0[ri.isz+isx*nz]-=fgetrick(it*ri.dt,ri.freq);

		step_forward(rp0,rp1,v,ri.nz,ri.nx,rz,rx,s);
		step_forward(p0,p1,v,ri.nz,ri.nx,rz,rx,s);

		/* swap pointers */
		tmp=rp0;
		rp0=rp1;
		rp1=tmp;
		tmp=p0;
		p0=p1;
		p1=tmp;
	    }

	    /* traces act as sources for backwards propagation 
	       of receiver field */
	    for (ix=0;ix<ntr;ix++) 
		rp1[ri.igz+(ri.igxbeg+ix)*nz]+=tr[ix*ri.nt+it];
		
	    /* option: read source field from movie file - seek is necessary
	       because the part of the file for each shot record is read
	       BACKWARDS */
	    if (NULL != ri.mfile) {
		sf_seek(ri.mfile,(nxz*it)*sizeof(float),SEEK_SET);
		sf_floatread(p0,nxz,ri.mfile);
	    }

	    /* imaging condition */
	    for (ix=1;ix<ri.nx-1;ix++) {
		for (iz=1;iz<nz-1;iz++) {
		    ioff=iz+ix*nz;
		    imag[ioff] += p0[ioff]*rp1[ioff];
		}
	    }
        
	    /* write receiver movie snap to file if necessary */
	    if (NULL != ri.rmfile) {
		sf_floatwrite(rp1,nxz,ri.rmfile);
		/* imf++; */
	    }
        
	    /* next t */
	}

	isx += ri.iskip;
	isrc++;
		
	/* next shot record */
    }

    /* write image file */
    sf_floatwrite(imag,nxz,ri.imfile);

    /* clean up */

    exit(0);
}     
