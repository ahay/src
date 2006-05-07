/* Diplet transform */
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
    int n1, n2, i3, n3, n12, n12p, np;
    bool inv;
    char *type;
    float *pp, *qq, eps, dp, p0;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if n, do inverse transform */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
   
    if (inv) {
	if (!sf_getint("np",&np)) sf_error("Need np=");
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");

	n3 = sf_leftsize(in,2);

	sf_putint(out,"n3",np);
	sf_putfloat(out,"o3",p0);
	sf_putfloat(out,"d3",dp);

	sf_putint(out,"n4",n3);
    } else {
	if (!sf_histint(in,"n3",&np)) sf_error("No n3= in input");
	if (!sf_histfloat(in,"o3",&p0)) sf_error("No o3= in input");
	if (!sf_histfloat(in,"d3",&dp)) sf_error("No d3= in input");

	n3 = sf_leftsize(in,3);

	sf_putint(out,"n3",n3);
    }

    n12 = n1*n2;
    n12p = n12*np;

    pp = sf_floatalloc(n12);
    qq = sf_floatalloc(n12p);
 
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* wavelet type (haar,linear) */

    diplet_init(n1,n2,np,p0,dp,inv,eps,type[0]);

    for (i3=0; i3 < n3; i3++) {
	if (inv) {	    
	    sf_floatread(pp,n12,in);
	} else {
	    sf_floatread(qq,n12p,in);
	}

	diplet_lop(inv,false,n12p,n12,qq,pp);

	if (inv) {
	    sf_floatwrite(qq,n12p,out);
	} else {
	    sf_floatwrite(pp,n12,out);
	}
    }
    
    exit(0);
}
