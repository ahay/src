/* Rice HPCSS forward modeling. */

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

#include "stepgeom.h"
#include "wavefun.h"

/* cfl number appropriate for scheme */
#define CFL 0.7


int main(int argc, char ** argv) {
    
    /* BEGIN DECLARATIONS */
    
    WINFO wi;        /* struct for command line input */
    
    /* workspace */
    
    float * v;       /* velocity field */
    float * p1;      /* pressure field, current time step */
    float * p0;      /* pressure field, last time step */
    
    float * tr;      /* storage for traces */
    float * tmp;     /* used to swap p1 and p0 */
    
    int ix, it;      /* counters */
    int isrc;        /* source counter */
    int imf;         /* movie frame counter */
    int isx;         /* source location, in units of dx */
    int nxz;         /* number of spatial grid points */
    int ntr;         /* number of traces */
    int nsam;        /* number of trace samples */
    int nsrc;        /* number of shots */
    float rz,rx,s;   /* precomputed coefficients */
    float vmax,vmin; /* max, min velocity values */
    int k_scale, num_cfl;

    /* geomB * arrB = (geomB*)0; */

    /* END DECLARATIONS */
    
    sf_init(argc,argv);
    
    /* read inputs from command line */
    getinputs(true,&wi);
    
    /* compute number of shots */
    nsrc = (wi.isxend-wi.isxbeg)/(wi.iskip); nsrc++;
    
    /* compute number of spatial grid points */
    nxz=wi.nx * wi.nz;
    
    /* compute number of traces, samples in each record */
    ntr=wi.igxend-wi.igxbeg+1;
    nsam=ntr*wi.nt;
    
    /* allocate, initialize p0, p1, v, traces */
    p0=sf_floatalloc(nxz);
    p1=sf_floatalloc(nxz);
    v =sf_floatalloc(nxz);
    tr=sf_floatalloc(nsam);
    
    /* read velocity */
    sf_floatread(v,nxz,wi.vfile);
    
    /* CFL, sanity checks */
    vmax=fgetmax(v,nxz);
    vmin=fgetmin(v,nxz);
    if (vmax*wi.dt>CFL*fmaxf(wi.dx,wi.dz)) {
	sf_warning("CFL criterion violated");
	sf_warning("vmax=%e dx=%e dz=%e dt=%e\n",vmax,wi.dx,wi.dz,wi.dt);
	sf_warning/*error*/("max permitted dt=%e\n",CFL*fmaxf(wi.dx,wi.dz)/vmax);
    }
    if (vmin<=0.0) 
	sf_error("min velocity nonpositive");
    
    /* only square of velocity array needed from here on */
    fsquare(v,nxz);
    
    /* precalculate some coefficients */
    rz=wi.dt*wi.dt/(wi.dz*wi.dz);
    rx=wi.dt*wi.dt/(wi.dx*wi.dx);
    s =2.0*(rz+rx);
    
    /* shot loop */
    isrc=0;
    isx=wi.isxbeg;

    /* INITIALIZE GEOMETRIC WAVE-CIRCLES 
    arrB = malloc(nxz * sizeof(geomB));
    construct_Boundaries(arrB, v, vmax, wi.dt,  wi.dz,  wi.dx, wi.nz, wi.nx); */

    k_scale = (int)(vmax * wi.dt / CFL / wi.dz) + 1;
    sf_warning("CFL k-scale = %d", k_scale);

    while (isx <= wi.isxend) {
	
	/* initialize pressure fields, traces */
	fzeros(p0,nxz);
	fzeros(p1,nxz);
	fzeros(tr,nsam);
	
	/* initialize movie frame counter */
	imf=0;
	
	/* time loop */
	sf_warning("t = 0: | p0 |_L1 = %g ",norm1(p0, wi.nz, wi.nx));

	for (it=0;it<wi.nt;it++) {
	    
	    /* construct next time step, overwrite on p0 */
	    
	    /* step_forward(p0,p1,v,wi.nz,wi.nx,rz,rx,s); 
	       step_forward_g(arrB, p0, p1, v, wi.nz, wi.nx, wi.dz, wi.dx, wi.dt, wi.isz, isx); */

	    num_cfl = step_forward_k(p0,p1,v,wi.nz,wi.nx,rz,rx,s, k_scale, wi.dt, CFL, wi.dx, wi.dz); 

	    sf_warning(" %g time |p0|_L1 num_cfl (per) = %g", it*100.f/wi.nt, (float)num_cfl*100.f/(float)(nxz));
	    /* norm1(p0, wi.nz, wi.nx)); */
	    
	    /* tack on source */
	    p0[wi.isz+isx*wi.nz]+=fgetrick(it*wi.dt,wi.freq);
	    sf_warning(" ... + getrick : | p0 |_L1 = %g ",norm1(p0, wi.nz, wi.nx));
	    
	    /* swap pointers */
	    tmp=p0;
	    p0=p1;
	    p1=tmp;
	    
	    /* store trace samples if necessary */
	    if (NULL != wi.tfile) 
		for (ix=0;ix<ntr;ix++) 
		    tr[ix*wi.nt+it]=p1[(wi.igxbeg+ix)*wi.nz+wi.igz];
	    
	    /* write movie snap to file if necessary */
	    if (NULL != wi.mfile && wi.nm && !(it%wi.nm)) {
		sf_floatwrite(p1,nxz,wi.mfile);
		imf++;
	    }
	    
	    /* next t */
	}

	/* write traces to file if necessary */
	if (NULL != wi.tfile) 
	    sf_floatwrite(tr,nsam,wi.tfile);
	
	isx += wi.iskip;
	isrc++;
    } 


    exit(0);
}
