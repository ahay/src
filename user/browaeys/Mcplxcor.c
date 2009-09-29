/* Statistical complex correlation for circular data. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <float.h>
#include <rsf.h>

#define SIGN(a) (a > 0 ? 1 : (a < 0 ? -1 : 0))

#include "cplxstat.h"

int main(int argc, char* argv[])
{
    int i,j,nt,nx;

    float c,p;
    float dt,dx,ot,ox;
    float **data,**cc,*t1,*t2; 

    sf_axis caxis,xaxis;
    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt))   sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&nx))   sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");

    /* output file parameters */
    caxis = sf_maxa(2,0,1);
    sf_oaxa(out,caxis,1);

    xaxis = sf_maxa (nx-1,ox,dx);
    sf_oaxa(out,xaxis,2);

    /* memory allocations */
    data = sf_floatalloc2(nt,nx);
    cc = sf_floatalloc2(2,nx-1);
    t1 = sf_floatalloc(nt);
    t2 = sf_floatalloc(nt);

    /* input in degree */
    sf_floatread(data[0],nt*nx,in);

    for (i = 0; i < nx-1; i++) { 

	for (j = 0; j < nt; j++) {

	    t1[j] = data[i][j];
	    t2[j] = data[i+1][j];
	
	}

	c = 0.0;  
	p = 0.0;

        /* cross-correlation modulus and phase */
	circ_corr(t1,t2,nt,&c,&p);

	cc[i][0] = p*180.0/SF_PI;            
	cc[i][1] = c;

    }

    sf_floatwrite(cc[0],2*(nx-1),out);

    exit (0);
}
