/* Local matching-radon operator. */
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

#include <math.h>

#include "matchoper.h"

int main(int argc, char* argv[])
{
    int m[SF_MAX_DIM], rect[SF_MAX_DIM];
    int mdim, nd, shift, niter;
    int i, i4, n4;
    float *data, *filt, anti, p1, pclip, *w=NULL;
    char key[6];
    bool verb;            
    bool par, freq, rho;  
    int nt;                   /* number of time samples */
    int nx;                   /* number of offsets */ 
    int np;                   /* number of slopes */
    float dt,t0;              /* time increment, starting time */
    float dp,p0;              /* slope increment, starting slope */
    float x0, dx, ox;         /* reference offset, increment, origin */
    sf_file in, out, weight=NULL;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* read input file parameters */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;
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

    if (!sf_getint("shift",&shift)) sf_error("Need shift=");
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

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

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getfloat("pclip",&pclip)) pclip=100.;

    mdim = sf_filedims(in,m);
    mdim = SF_MIN(mdim,2);
    data = sf_floatalloc(nt*nx);
    filt = sf_floatalloc(nt*nx*shift);

    if (NULL != sf_getstring ("weight")) {
	weight = sf_output("weight");
	w = sf_floatalloc(nt*np);
	sf_putint(weight,"n2",np);
	sf_putfloat(weight,"d2",dp);
	sf_putfloat(weight,"o2",p0);
	sf_putstring(weight,"label2","Weight");
	sf_putstring(weight,"unit2","");
    } else {
	weight = NULL;
	w = NULL;
    }

    sf_shiftdim(in, out, 2);
    sf_putint(out,"n3",shift);
    sf_putfloat(out,"d3",1);
    sf_putfloat(out,"o3",0);
    sf_putstring(out,"label3","Filter");
    sf_putstring(out,"unit3","");

    nd = 1;
    for (i=0; i < mdim; i++) {
	nd *= m[i];
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
    }

    if (nt*nx != nd) sf_error("Wrong data dimensions!");

    n4 = sf_leftsize(in,2);

    for (i=mdim; i < SF_MAX_DIM; i++) {
	m[i] = 1;
    }
    
    matchoper_init (shift,mdim,niter,m,rect,verb,nt,nx,np,dt,
		    t0,dx,ox,x0,dp,p0,par,freq,rho,anti,p1,pclip);

    for (i4 = 0; i4 < n4; i4++) { /* loop over CMPs */
	if(verb) sf_warning("i=%d of %d",i4+1,n4);

	sf_floatread(data,nt*nx,in);

	matchoper_lop(filt,data,w);

	sf_floatwrite(filt,nt*nx*shift,out);

	if (NULL!=weight) {
	    sf_floatwrite(w,nt*np,weight);
	}
    }

    exit(0);
}

/* 	$Id$	 */
