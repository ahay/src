/* Time-space-domain Radon transform (slant stack) 

April 2015 program of the month:
http://ahay.org/blog/2015/04/21/program-of-the-month-sfslant/
*/
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

#include <rsf.h>
#include "slant.h"

int main(int argc, char* argv[])
{
    int nt, n3, nx, np, i3, ntx, ntp;
    bool adj,verb,rho;
    float o1,d1, x0,dx, anti, p0,dp,p1;
    float *cmp=NULL, *vscan=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);

    in  = sf_input ("in");
    out = sf_output("out");

    if (!sf_histint  (in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    n3 = sf_leftsize(in,2);

    if (!sf_getbool ("verb",&verb)) verb=false; /* verbosity flag */
    if (!sf_getbool ( "adj",&adj )) adj=false;  /* adjoint flag */
    if (!sf_getbool ( "rho",&rho )) rho=true;   /* rho filtering */
    if (!sf_getfloat("anti",&anti)) anti=1.;    /* antialiasing */

    if (adj) {
	if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histint  (in,"n2",&nx)) sf_error("No n2= in input");

	/* specify slope axis */
	if (!sf_getint("np",&np)) sf_error("Need np=");
	/* number of p values (if adj=y) */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
	/* p sampling (if adj=y) */
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
	/* p origin (if adj=y) */

	sf_putint  (out,"n2",np);
	sf_putfloat(out,"d2",dp);
	sf_putfloat(out,"o2",p0);
    } else { /* modeling */
	if (!sf_histint  (in,"n2",&np)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");

	if (!sf_getfloat("x0",&x0)) sf_error("Need x0=");
	/* offset origin */
	if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
	/* offset sampling */
	if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
	/* number of offsets */

	sf_putint  (out,"n2",nx);
	sf_putfloat(out,"o2",x0);
	sf_putfloat(out,"d2",dx);
    }

    if (!sf_getfloat("p1",&p1)) p1=0.;
    /* reference slope */

    ntx = nt*nx;
    ntp = nt*np;

    cmp   = sf_floatalloc(ntx);
    vscan = sf_floatalloc(ntp);

    slant_init (true, rho, x0, dx, nx, p0, dp, np, o1, d1, nt, p1, anti);

    for (i3=0; i3 < n3; i3++) { 
	if(verb) sf_warning("i=%d of %d",i3+1,n3);
	if( adj) {
	    sf_floatread(  cmp,ntx,in);
	} else {
	    sf_floatread(vscan,ntp,in);
	}

	slant_lop(adj,false,ntp,ntx,vscan,cmp);

	if( adj) {
	    sf_floatwrite(vscan,ntp,out);
	} else {
	    sf_floatwrite(  cmp,ntx,out);
	} 
    }

    exit(0);
}
