/* 2-D (inline) trace interpolation to a denser grid using PWD.

It may be necessary to bandpass the data before and after dealiasing 
to ensure that the temporal spectrum is banded. Rule of thumb: if 
max(jx,jy)=N, the temporal bandwidth should be 1/N of Nyquist.
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

#include <math.h>

#include <rsf.h>

#include "predict.h"

int main (int argc, char *argv[])
{
    bool both;
    int n1, n2, n3, n12, m2, m12, i1, i2, i3, order;
    float d2, eps;
    float **u1, **p, **uu, **p2;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if(SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3=sf_leftsize(in,2);

    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    m2 = n2*2-1; 
    if (n2 > 1) d2 /= 2;

    n12 = n1*n2;
    m12 = n1*m2;
    
    sf_putint(out,"n2",m2);
    sf_putfloat(out,"d2",d2);

    u1 = sf_floatalloc2(n1,n2);
    uu = sf_floatalloc2(n1,m2);

    if (!sf_getbool("both",&both)) both=false;
    /* if use left and right slopes */

    p = sf_floatalloc2(n1,n2);
    p2 = both? sf_floatalloc2(n1,n2): NULL;

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init (n1,1,eps*eps,order,0,true);

    for (i3=0; i3 < n3; i3++) {
	/* get data */
	sf_floatread(u1[0],n12,in);

	/* get slopes */
	sf_floatread(p[0],n12,dip);
	if (both) sf_floatread(p2[0],n12,dip);

	/* interpolate inline */
	for (i2=0; i2 < n2-1; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		uu[2*i2][i1] = u1[i2][i1];
	    }
	    predict2_step(true,both,u1[i2],u1[i2+1],
			  p[i2],both? p2[i2+1]:p[i2+1],
			  uu[2*i2+1]);
	}
	for (i1=0; i1 < n1; i1++) {
	    uu[2*n2-2][i1] = u1[n2-1][i1];
	}
    
	sf_floatwrite(uu[0],m12,out);
    }

    exit (0);
}

/* 	$Id: Mdealias2.c 1713 2006-03-03 08:21:29Z fomels $	 */
