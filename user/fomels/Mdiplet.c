/* 2-D Seislet frame */
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

#include "diplet.h"

int main(int argc, char *argv[])
{
    int i, n1, n2, i3, n3, n12, np, n12p, ncycle, niter;
    bool inv, verb;
    char *type;
    float *pp, *qq, ***dd, eps, perc, scale;
    sf_file in, out, dip;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dips");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    if (!sf_histint(dip,"n3",&np)) np=1;
    n12p = n12*np;
    scale = 1./np;

    pp = sf_floatalloc(n12);
    qq = sf_floatalloc(n12p);
    dd = sf_floatalloc3(n1,n2,np);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getint("ncycle",&ncycle)) ncycle=0;
    /* number of iterations */

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of Bregman iterations */

    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    if (inv) {
	n3 = sf_leftsize(in,3);
	sf_putint(out,"n3",n3);
    } else {
	n3 = sf_leftsize(in,2);
	sf_putint(out,"n3",np);
	sf_putint(out,"n4",n3);
    } 

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* wavelet type (haar,linear,biorthogonal), default is linear */
    
    diplet_init(n1,n2,np,dd,true,eps,type[0]);
    
    for (i3=0; i3 < n3; i3++) {
	sf_floatread(dd[0][0],n12p,dip);	
	
	if (!inv) {
    	    sf_floatread(pp,n12,in);
	    diplet_lop(true,false,n12p,n12,qq,pp);
	    if (ncycle > 0) 
		sf_sharpinv(diplet_lop,
			    scale,niter,ncycle,perc,verb,n12p,n12,qq,pp); 
	    sf_floatwrite(qq,n12p,out);

	} else {
	    sf_floatread(qq,n12p,in);
	    diplet_lop(false,false,n12p,n12,qq,pp);

	    for (i=0; i < n12; i++) {
		pp[i] *= scale;
	    }
	    sf_floatwrite(pp,n12,out);
	}
    }
    
    exit(0);
}
