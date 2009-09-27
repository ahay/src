/* 2-D trace interpolation to a denser grid using PWD.

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
#include "interpd.h"

int main (int argc, char *argv[])
{
    int n1, n2, n3, n12, m2, m3, m12, i1, i2, i3;
    float d2, d3, ***u1, ***uu1, ***u2, ***uu2, ***p, ***q, ***q2, ***uu;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if(SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) n3=1;

    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;

    m2 = n2*2-1; 
    if (n2 > 1) d2 /= 2;

    m3 = n3*2-1; 
    if (n3 > 1) d3 /= 2;

    n12 = n1*n2*n3;
    m12 = n1*m2*m3;
    
    sf_putint(out,"n2",m2);
    sf_putint(out,"n3",m3);
 
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"d3",d3);

    u1 = sf_floatalloc3(n1,n2,n3);
    uu = sf_floatalloc3(n1,m2,m3);

    p = sf_floatalloc3(n1,n2,n3);

    q = sf_floatalloc3(n1,n2,n3);
    q2 = sf_floatalloc3(n1,n3,m2);
    
    uu1 = sf_floatalloc3(n1,m2,n3);

    u2 = (float***) sf_alloc(m2,sizeof(float**));
    uu2 = (float***) sf_alloc(m2,sizeof(float**));

    for (i2=0; i2 < m2; i2++) {
	u2[i2]  = (float**) sf_alloc(m3,sizeof(float*));
	for (i3=0; i3 < m3; i3++) {
	    u2[i2][i3] = uu[i3][i2];
	} 
	uu2[i2] = (float**) sf_alloc(n3,sizeof(float*));
	for (i3=0; i3 < n3; i3++) {
	    uu2[i2][i3] = uu1[i3][i2];
	}
    }

    interp_init (n1, 0.0001);

    /* get data */
    sf_floatread(u1[0][0],n12,in);

    /* get slopes */
    sf_floatread(p[0][0],n12,dip);
    if (n3 > 1) sf_floatread(q[0][0],n12,dip);

    /* interpolate inline */
    for (i3=0; i3 < n3; i3++) {
	interp2(n2,u1[i3],uu1[i3],p[i3]);
    }
    
    /* extend crossline slope */    
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2-1; i2++) {
	    q2[2*i2][i3] = q[i3][i2];
	    for (i1=0; i1 < n1; i1++) {
		q2[2*i2+1][i3][i1] = 0.5*(q[i3][i2][i1]+q[i3][i2+1][i1]);
	    }
	}
    }

    /* interpolate crossline */
    for (i2=0; i2 < m2; i2++) {
	interp2(n3,uu2[i2],u2[i2],q2[i2]);
    }
    
    sf_floatwrite(uu[0][0],m12,out);
    sf_close();
    exit (0);
}

/* 	$Id: Mdealias2.c 1713 2006-03-03 08:21:29Z fomels $	 */
