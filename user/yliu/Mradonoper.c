/* Linear Radon operator. */
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
#include "radonoper.h"       

int main (int argc, char **argv)
{
    bool verb;                /* verbosity flag */
    int nt;                   /* number of time samples */
    int nx;                   /* number of offsets */ 
    int np;                   /* number of slopes */
    int nc,ic;                /* number of CMP gathers, CMP gather counter */
    bool adj, par;            /* adjoint, parabolic */
    float dt,t0;              /* time increment, starting time */
    float dp,p0;              /* slope increment, starting slope */
    float x0, dx, ox;         /* reference offset, increment, origin */
    float *x, *y;             /* input and output */
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("adj",&adj)) adj=true;
    /* if y, perform adjoint operation */

    if (!sf_getbool ("verb",&verb)) verb=false; /* verbosity flag */

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;

    if (adj) { 
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

	/* specify slope axis */

	if (!sf_getint  ("np",&np)) sf_error("Need np=");
	/* number of p values (if adj=y) */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
	/* p sampling (if adj=y) */
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
        /* p origin (if adj=y) */

	sf_putint(  out,"n2",np);
	sf_putfloat(out,"d2",dp);
	sf_putfloat(out,"o2",p0);

    } else { /* modeling */

	if (!sf_histint  (in,"n2",&np)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");

	if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
	/* number of offsets (if adj=n) */
	
	sf_putint(out,"n2",nx);
    }

    nc = sf_leftsize(in,2);

    if (adj) {
	if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    } else {
	if (!sf_getfloat("ox",&ox)) sf_error("Need ox=");
	if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
    }
    
    if (!adj) {
	sf_putfloat(out,"o2",ox);
	sf_putfloat(out,"d2",dx);
    }
    
    if (!sf_getbool("parab",&par)) par=false;
    /* if y, parabolic Radon transform */
    if (!sf_getfloat("x0",&x0)) x0=1.;
    /* reference offset */

    x = sf_floatalloc (nt*np);
    y = sf_floatalloc (nt*nx);

    radonoper_init(nt,dt,t0,nx,dx,ox,x0,np,dp,p0,par);

    for (ic = 0; ic < nc; ic++) { /* loop over CMPs */
	if(verb) sf_warning("i=%d of %d",ic+1,nc);
	if (adj) {
	    sf_floatread(y,nt*nx,in);
	} else { /* modeling */
	    sf_floatread(x,nt*np,in);
	}

	radonoper_lop (adj, false, nt*np, nt*nx, x, y);

	if (adj) {
	    sf_floatwrite(x,nt*np,out);
	} else { /* modeling */
	    sf_floatwrite(y,nt*nx,out);
       	}
    } /* loop over CMPs */

    radonoper_close();

    exit (0);
}

/* 	$Id$	 */
