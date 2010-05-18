/* Test linear cascading matching-Radon operator. */
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
#include "nmrt.h"

int main (int argc, char **argv)
{
    bool verb;                /* verbosity flag */
    int nt;                   /* number of time samples */
    int nx;                   /* number of offsets */ 
    int np;                   /* number of slopes */
    int nc,ic;                /* number of CMP gathers, CMP gather counter */
    bool adj, inv, par, freq; /* adjoint, inverse, parabolic, domain flag */
    int nf1, nf2, nf;         /* number of filter size */
    int niter;                /* number of conjugate gradient iteration */
    float dt,t0;              /* time increment, starting time */
    float dp,p0;              /* slope increment, starting slope */
    float x0, dx, ox;         /* reference offset, increment, origin */
    float *x, *y, *filt;      /* input and output */
    bool rho;
    float anti, p1;
    sf_file in, out, fil;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    fil = sf_input("filt");

    if (!sf_getbool("adj",&adj)) adj=false;
    /* if y, perform adjoint operation */

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, perform inverse operation */

    if (!sf_getbool ("verb",&verb)) verb=false; /* verbosity flag */

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

    if (!sf_histint(fil,"n1",&nf1)) sf_error("No nf1= in filt");
    if (!sf_histint(fil,"n2",&nf2)) sf_error("No nf2= in filt");
    if (!sf_histint(fil,"n3",&nf)) sf_error("No nf3= in filt");

    if (nt != nf1 || nx != nf2 ) sf_error("Nee n1==nf1 && n2==nf2");

    /* specify slope axis */
    
    if (!sf_getint  ("np",&np)) sf_error("Need np=");
    /* number of p values */
    if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
    /* p sampling */
    if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
    /* p origin */

    nc = sf_leftsize(in,2);

    if (!sf_getint("niter",&niter)) niter=100;;
    /* number of conjugate gradient iteration */
   

    if (!sf_getbool("freq",&freq)) freq=true;
    /* if y, parabolic Radon transform */

    if (!sf_getbool("parab",&par)) par=false;
    /* if y, parabolic Radon transform, only when freq=y */
    if (!sf_getfloat("x0",&x0)) x0=1.;
    /* reference offset */

    if (!sf_getbool ( "rho",&rho )) rho=true;
    /* rho filtering, only when freq=n */
    if (!sf_getfloat("anti",&anti)) anti=1.;
    /* antialiasing, only when freq=n */
    if (!sf_getfloat("p1",&p1)) p1=0.;
    /* reference slope, only when freq=n */
    
    x = sf_floatalloc (nt*nx);
    y = sf_floatalloc (nt*nx);
    filt = sf_floatalloc (nt*nx*nf);

    for (ic = 0; ic < nc; ic++) { /* loop over CMPs */
	if(verb) sf_warning("i=%d of %d",ic+1,nc);

	sf_floatread(filt,nt*nx*nf,fil);
	nmrt_init (nt,nx,nf,np,dt,t0,dx,ox,x0,dp,p0,par,freq,rho,anti,p1,filt);
	if (!inv) {
	    if (adj) {
		sf_floatread(y,nt*nx,in);
	    } else { /* modeling */
		sf_floatread(x,nt*nx,in);
	    }
	    
	    nmrt_oper (adj, false, nt*nx, nt*nx, x, y);
	    if (adj) {
		sf_floatwrite(x,nt*nx,out);
	    } else { /* modeling */
		sf_floatwrite(y,nt*nx,out);
	    }
	} else {
	    sf_floatread(y,nt*nx,in);
	    nmrt_lop(niter,x,y,verb);
	    sf_floatwrite(x,nt*nx,out);
	}

    } /* loop over CMPs */

    exit (0);
}

/* 	$Id: Mtestcasoper.c 5959 2010-05-14 04:11:09Z yang_liu $	 */
