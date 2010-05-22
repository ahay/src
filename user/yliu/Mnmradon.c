/* Nonstationary-matching Radon transform. */
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
#include "radonoper.h"
#include "slant.h"

int main (int argc, char **argv)
{
    bool verb;                /* verbosity flag */
    int nt;                   /* number of time samples */
    int nx;                   /* number of offsets */ 
    int np;                   /* number of slopes */
    int nw1, nw2;
    int nc,ic,i;              /* number of CMP gathers, CMP gather counter */
    bool adj, inv, par, freq; /* adjoint, inverse, parabolic, domain flag */
    int nf1, nf2, nf;         /* number of filter size */
    float dt,t0;              /* time increment, starting time */
    float dp,p0;              /* slope increment, starting slope */
    float x0, dx, ox;         /* reference offset, increment, origin */
    float *x, *y, *filt;      /* input and output */
    bool rho;
    float anti, p1, *w=NULL;
    sf_file in, out, fil, weight=NULL;

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

    if (!sf_histint(fil,"n1",&nf1)) sf_error("No nf1= in filt");
    if (!sf_histint(fil,"n2",&nf2)) sf_error("No nf2= in filt");
    if (!sf_histint(fil,"n3",&nf)) sf_error("No nf3= in filt");


    if (!adj) {
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
	
	/* specify slope axis */
	if (!sf_getint  ("np",&np)) sf_error("Need np=");
	/* number of p values */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
	/* p sampling */
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
	/* p origin */
    } else {
	if (!sf_histint(in,"n2",&np)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
	
	/* specify space axis */
	if (!sf_getint  ("nx",&nx)) sf_error("Need nx=");
	/* number of x values */
	if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
	/* x sampling */
	if (!sf_getfloat("ox",&ox)) sf_error("Need ox=");
	/* x origin */
    }

    if (nf1!=nt || nf2!= nx) sf_error("Need nf1==nt && nf2==nx");

    nc = sf_leftsize(in,2);

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
    
    x = sf_floatalloc (nt*np);
    y = sf_floatalloc (nt*nx);
    filt = sf_floatalloc (nt*nx*nf);

    if (NULL != sf_getstring ("weight")) {
	weight = sf_input("weight");
	w = sf_floatalloc(nt*np);
	if (!sf_histint(in,"n1",&nw1)) sf_error("No n1= in weight");
	if (!sf_histint(in,"n2",&nw2)) sf_error("No n2= in weight");
	if (nw1 != nt || nw2 != np) sf_error("Need nw1==nt && nw2==np");
    } else {
	weight = NULL;
	w = NULL;
    }


    for (ic = 0; ic < nc; ic++) { /* loop over CMPs */
	if(verb) sf_warning("i=%d of %d",ic+1,nc);

	sf_floatread(filt,nt*nx*nf,fil);

	if (!adj) {
	    sf_floatread(y,nt*nx,in);
	    nmrt_init (nt,nx,nf,np,dt,t0,dx,ox,x0,dp,p0,
		       par,freq,rho,anti,p1,filt);

	    nmrt_oper (adj, false, nt*np, nt*nx, x, y);

	    if (NULL != weight) {
		sf_floatread(w,nt*np,weight);
		for (i=0; i < nt*np; i++) {
		    x[i] *= w[i];
		}
	    }

	    sf_putint(out,"n2",np);
	    sf_putfloat(out,"d2",dp);
	    sf_putfloat(out,"o2",p0);
	    sf_putstring(out,"label2","P");
	    sf_putstring(out,"unit2","");
	    	    
	    sf_floatwrite(x,nt*np,out);

	} else {
	    sf_floatread(x,nt*np,in);

	    if (inv) {
		if (freq) {
		    radonoper_init (nt,dt,t0,nx,dx,ox,
				    x0,np,dp,p0,par);
		    radonoper_lop (false, false, nt*np, nt*nx, x, y);
		    
		} else {
		    slant_init (true,rho,x0,dx,nx,p0,dp,
				np,t0,dt,nt,p1,anti);
		    slant_lop(false,false,nt*np,nt*nx,x,y);
		}
	    } else {
		nmrt_init (nt,nx,nf,np,dt,t0,dx,ox,x0,dp,p0,
			   par,freq,rho,anti,p1,filt);
	
		nmrt_oper (adj, false, nt*np, nt*nx, x, y);
	    }
	    sf_putint(out,"n2",nx);
	    sf_putfloat(out,"d2",dx);
	    sf_putfloat(out,"o2",ox);
	    sf_putstring(out,"label2","Distance");
	    sf_putstring(out,"unit2","");

	    sf_floatwrite(y,nt*nx,out);
	}

    } /* loop over CMPs */

    exit (0);
}

/* 	$Id$	 */
