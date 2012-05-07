/* 2-D omnidirectional plane-wave destruction */
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

static void lagrange(float p, float *f)
{
    f[0]=0.5*p*(p-1.0);
    f[1]=1.0-p*p;
    f[2]=0.5*p*(p+1.0);
}

int main(int argc, char *argv[])
{
    int i1, i2, n1, n2, n12;
    float **dat, **ang, **res, filt[3], *trace;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    dip = sf_input("dip");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    dat = sf_floatalloc2(n1,n2);
    ang = sf_floatalloc2(n1,n2);
    res = sf_floatalloc2(n2,n1);

    sf_floatread(dat[0],n12,inp);
    sf_floatread(ang[0],n12,dip);

    for (i2=0; i2 < n2; i2++) {
	res[0][i2] = res[n1-1][i2] = 0.;
	trace = dat[i2];
	for (i1=1; i1 < n1-1; i1++) {
	    lagrange(sinf(ang[i2][i1]),filt);

	    res[i1][i2] = filt[0]*trace[i1-1]+filt[1]*trace[i1]+filt[2]*trace[i1+1];
	}
	dat[i2][0] = dat[i2][n1-1] = 0.;
    }
    for (i1=0; i1 < n1; i1++) {
	dat[0][i1] = dat[n2-1][i1] = 0.;
	trace = res[i1];
	for (i2=1; i2 < n2-1; i2++) {
	    lagrange(cosf(ang[i2][i1]),filt);

	    dat[i2][i1] -= filt[0]*trace[i1-1]+filt[1]*trace[i1]+filt[2]*trace[i1+1];
	}
    }

    sf_floatwrite(dat[0],n12,out);

    exit(0);
}
