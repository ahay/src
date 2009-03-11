/* Complex coherency based on circular statistics cross-correlation. */
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

#include <rsf.h>
#include <math.h>

#define SIGN(a) (a > 0 ? 1 : (a < 0 ? -1 : 0))

void circ_mean (float *d, int n, float *v, float *t)
{
    int i;
    float r,c,s;

    r = 1.0;

    c = 0.0; s = 0.0;

    for (i = 0; i < n; i++) {
	c += cos(r*d[i]);
	s += sin(r*d[i]);
    }

    c /= n; s /= n;

    *v = 1.0 - (c*c + s*s);

    *t = atan2(s,c);
    if (*t < 0.0) *t += 2.0*SF_PI;

    return;
}

void circ_corr (float *d1, float *d2, int n, float *m, float *p)
{
    int i;

    float r,v1,v2,t1,t2;
    float c,s,d,r1r2,rm,cm;

    r = 1.0;

    circ_mean(d1,n,&v1,&t1);
    circ_mean(d2,n,&v2,&t2);

    c = 0.0; s = 0.0;

    for (i = 0; i < n; i++) {
	d = r*(d1[i] - d2[i]);
	c += cos(d);
	s += sin(d);
    }

    c /= n; s /= n;

    r1r2 = sqrt(1.0 - v1)*sqrt(1.0 - v2);

    rm = c - r1r2*cos(t1 - t2);
    cm = s - r1r2*sin(t1 - t2);

    /* modulus */
    *m = sqrt((rm*rm + cm*cm)/(v1*v2));

    /* phase */
    *p = atan2(cm,rm);
    if (*p < 0.0) *p += 2.0*SF_PI;

    return;
}


int main(int argc, char* argv[])
{
    int i,k,j;
    int nt,nx,nw,ntcc;

    float dt,dx,ot,ox;
    float t,c,p;

    float *t1,*t2;
    float **data,**cc,**cd; 

    sf_axis xaxis,taxis;
    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");

    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");

    if (!sf_getint("nw",&nw)) sf_error("No nw= in input");
    /* half time-window size */

    ntcc = nt - 2*nw;

    /* output file parameters */
    xaxis = sf_maxa (nx-1, ox, dx);
    sf_oaxa(out,xaxis,1);

    taxis = sf_maxa(ntcc,nw,dt);
    sf_oaxa(out,taxis,2);

    data = sf_floatalloc2(nt,nx);
    cc = sf_floatalloc2(ntcc,nx-1);
    cd = sf_floatalloc2(ntcc,nx-1);

    sf_floatread(data[0],nt*nx,in);

    t1 = sf_floatalloc(2*nw+1);
    t2 = sf_floatalloc(2*nw+1);

    for (i = 0; i < 2*nw+1; i++) {
	t1[i] = 0.0;
	t2[i] = 0.0;
    }

    for (i = 0; i < nx-1; i++) { 

	for (k = nw; k < nt-nw; k++) { 

	    t = ot + k*dt; /* time location */

	    for (j = -nw; j < nw+1; j++) {
                /* time window */
		t1[j] = data[i][k+j];
		t2[j] = data[i+1][k+j];
	    }

	    c = 0.0;  /* coherency = modulus*/
	    p = 0.0;  /* dip = phase */

            /* cross-correlation modulus and phase */
	    circ_corr(t1,t2,2*nw+1,&c,&p);

            /* coherency = modulus*/
	    cc[i][k] = c;

            /* dip = phase */
	    cd[i][k] = p;

	}

    }

    sf_floatwrite (cc[0],ntcc*(nx-1),out);
    
    exit (0);
}


