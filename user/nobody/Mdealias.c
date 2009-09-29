/* Trace interpolation to a denser XY grid using PWD.

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

#include <rsf.h>

#include "dip3.h"
#include "interpd.h"

int main (int argc, char *argv[])
{
    int n1, n2, n3, n12, niter, nf, i, i2, i1, m1, m2, m3, m12, nj1, nj2;
    int rect[3], liter;
    bool verb, sign;
    float d2, d3, p0, q0, *t, *tt;
    float ***u1, ***uu1, ***p, ***q, ***uu, ***u2, ***uu2, ***q2;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) n3=1;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d3",&d3)) d3=d2;

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of dip-estimation iterations */
    if (!sf_getint("liter",&liter)) liter=100;
    /* number of linear iterations */
    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    if (!sf_getint("rect3",&rect[2])) rect[2]=1;
    /* dip smoothness */

    if (!sf_getint("order",&nf)) nf=1;
    /* [1,2,3] dip filter accuracy */
    if (nf < 1 || nf > 3) sf_error ("accuracy must be between 1 and 3");

    if(SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if(!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if(!sf_getbool("sign",&sign)) sign = false;
    /* if y, keep dip sign constant */

    if (!sf_getfloat("p0",&p0)) p0=0.;
    /* initial in-line dip */
    if (!sf_getfloat("q0",&q0)) q0=0.;
    /* initial cross-line dip */

    if (!sf_getint("nj1",&nj1)) nj1=2;
    /* in-line antialiasing */
    if (!sf_getint("nj2",&nj2)) nj2=2;
    /* cross-line antialiasing */

    m1 = n1;
    m2 = n2*2-1; if (n2 > 1) d2 /= 2;
    m3 = n3*2-1; if (n3 > 1) d3 /= 2;
  
    n12 = n1*n2*n3;
    m12 = m1*m2*m3;

    sf_putint(out,"n2",m2); 
    sf_putfloat(out,"d2",d2);
    sf_putint(out,"n3",m3); 
    sf_putfloat(out,"d3",d3);

    u1 = sf_floatalloc3(n1,n2,n3);
    t = sf_floatalloc(n12);

    tt = sf_floatalloc(m12);
    uu = sf_floatalloc3(m1,m2,m3);
    u2 = sf_floatalloc3(m1,m3,m2);

    p = sf_floatalloc3(n1,n2,n3);
    q = sf_floatalloc3(n1,n2,n3);

    uu1 = sf_floatalloc3(n1,m2,n3);
    uu2 = sf_floatalloc3(n1,n3,m2);
    q2  = sf_floatalloc3(n1,n3,m2);

    free (uu2[0][0]);
    for (i1=0; i1 < m2; i1++) {
	for (i2=0; i2 < n3; i2++) {
	    uu2[i1][i2] = uu1[i2][i1];
	}
    }
 
    free (u2[0][0]);
    for (i1=0; i1 < m2; i1++) {
	for (i2=0; i2 < m3; i2++) {
	    u2[i1][i2] = uu[i2][i1];
	}
    }
 
    dip3_init(n1, n2, n3, rect, liter, sign);
    interp_init (n1, 0.0001);

    sf_floatread(u1[0][0],n12,in);

    if (verb) sf_warning("Estimating slopes...");
    
    for (i2=0; i2 < n3; i2++) {
	for (i1=0; i1 < n2; i1++) {
	    for (i=0; i < n1; i++) {
		p[i2][i1][i]=p0;
		q[i2][i1][i]=q0;
	    }
	}
    }

    dip3(1, niter, nf, nj1, verb, u1, p, NULL);
    dip3(2, niter, nf, nj2, verb, u1, q, NULL);

    if (verb) sf_warning("Expanding slopes...");

/*
    if (nj1 > 1) {
	for (i=0; i < n12; i++) {
	    p[0][0][i] /= nj1;
	}
    }

    if (nj2 > 1) {
	for (i=0; i < n12; i++) {
	    q[0][0][i] /= nj2;
	}
    }
*/ 

    for (i2=0; i2 < n3; i2++) {
	for (i1=0; i1 < n2-1; i1++) {
	    q2[2*i1][i2] = q[i2][i1];
	    for (i=0; i < n1; i++) {
		q2[2*i1+1][i2][i] = 0.5*(q[i2][i1][i]+q[i2][i1+1][i]);
	    }
	}
    }


    if (verb) sf_warning("Interpolating...");

    for (i2=0; i2 < n3; i2++) {
	interp2(n2,u1[i2],uu1[i2],p[i2]);
    }
    
    for (i1=0; i1 < m2; i1++) {
	interp2(n3,uu2[i1],u2[i1],q2[i1]);
    }
    
    sf_floatwrite(uu[0][0],m12,out);

    exit (0);
}

/* 	$Id$	 */
