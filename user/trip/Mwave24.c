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

#include "step24.h"
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
    /* int nz;          local number of gridpoints */
    int ntr;         /* number of traces */
    int nsam;        /* number of trace samples */
    int nsrc;        /* number of shots */
    float rz,rx,s;   /* precomputed coefficients */
    float vmax,vmin; /* max, min velocity values */
    /* float two;       two */
    
    /* END DECLARATIONS */
    
    sf_init(argc,argv);
    
    /* read inputs from command line */
    getinputs(true,&wi);
    
    /* compute number of shots */
    nsrc = (wi.isxend-wi.isxbeg)/(wi.iskip); nsrc++;
    
    /* compute number of spatial grid points */
    /* 2,4 FEATURE: nz->nz+2, nx->nx+2 */
    nxz=(wi.nx+2) * (wi.nz+2);
    
    /* compute number of traces, samples in each record */
    ntr=wi.igxend-wi.igxbeg+1;
    nsam=ntr*wi.nt;
    
    /* allocate, initialize p0, p1, v, traces */
    p0=sf_floatalloc(nxz);
    p1=sf_floatalloc(nxz);
    v =sf_floatalloc(nxz);
    tr=sf_floatalloc(nsam);
    
    /* read velocity */
    /* 2,4 FEATURE: nz->nz+2, etc: read interior from file, extend by const */
    /* read into interior in x */
    for (ix=1;ix<wi.nx+1;ix++) {
	/* read into interior of each column */
	sf_floatread(v+ix*(wi.nz+2)+1,wi.nz,wi.vfile);
	/* extend by const */
	v[ix*(wi.nz+2)]=v[ix*(wi.nz+2)+1];
	v[(ix+1)*(wi.nz+2)-1]=v[(ix+1)*(wi.nz+2)-2];
    }
    /* extend to 1st and last cols */
    for (ix=0;ix<wi.nz+2;ix++) {
	v[ix]=v[ix+wi.nz+2];
	v[ix+(wi.nx+1)*(wi.nz+2)]=v[ix+(wi.nx)*(wi.nz+2)];
    }
    
    /* CFL, sanity checks */
    vmax=fgetmax(v,nxz);
    vmin=fgetmin(v,nxz);
    if (vmax*wi.dt>CFL*fmaxf(wi.dx,wi.dz)) {
	sf_warning("CFL criterion violated");
	sf_warning("vmax=%e dx=%e dz=%e dt=%e\n",vmax,wi.dx,wi.dz,wi.dt);
	sf_error("max permitted dt=%e\n",CFL*fmaxf(wi.dx,wi.dz)/vmax);
    }
    if (vmin<=0.0) 
	sf_error("min velocity nonpositive");
    
    /* only square of velocity array needed from here on */
    fsquare(v,nxz);
    
    /* precalculate some coefficients */
    rz=wi.dt*wi.dt/(wi.dz*wi.dz);
    rx=wi.dt*wi.dt/(wi.dx*wi.dx);
    s =2.0*(rz+rx);
/*    two=2.0;  */
/*    nz=wi.nz; */
    
    /* shot loop */
    isrc=0;
    isx=wi.isxbeg;
    while (isx <= wi.isxend) {
	
	/* initialize pressure fields, traces */
	fzeros(p0,nxz);
	fzeros(p1,nxz);
	fzeros(tr,nsam);
	
	/* initialize movie frame counter */
	imf=0;
	
	/* time loop */
	for (it=0;it<wi.nt;it++) {
	    
	    /* construct next time step, overwrite on p0 */
	    
	    /* 2,4 FEATURE: nz->nz+2, nx->nx+2 */
	    step24_forward(p0,p1,v,wi.nz+2,wi.nx+2,rz,rx,s);
	    
	    /* tack on source */
	    /* 2,4 FEATURE: nz->nz+2, isz->isz+1, isx->isx+1 */
	    p0[(wi.isz+1)+(isx+1)*(wi.nz+2)]+=fgetrick(it*wi.dt,wi.freq);
	    
	    /* swap pointers */
	    tmp=p0;
	    p0=p1;
	    p1=tmp;
	    
	    /* store trace samples if necessary */
	    if (NULL != wi.tfile) 
		for (ix=0;ix<ntr;ix++) 
                    /* 2,4 FEATURE: nz->nz+2, igz->igz+1 */
		    tr[ix*wi.nt+it]=p1[(wi.igxbeg+ix+1)*(wi.nz+2)+wi.igz+1];
	    
	    /* write movie snap to file if necessary */
	    if (NULL != wi.mfile && wi.nm && !(it%wi.nm)) {
                /* 2,4 FEATURE: write interior of each interior column */
		for (ix=1;ix<wi.nx+1;ix++) 
		    sf_floatwrite(p1+ix*(wi.nz+2)+1,wi.nz,wi.mfile);
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
