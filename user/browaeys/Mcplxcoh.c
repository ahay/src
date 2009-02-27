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

    r = SF_PI/180.0;

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

    r = SF_PI/180.0;

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
    int i,j,k,it;

    int nt,nx,itype,nts,ntw,ntcc;

    float dt,dx,ot,ox;
    float t,ts,cm,pm,m,phi;
    float *t1,*t2;
    float **trc,**cc; 

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
    
    if (!sf_getint("itype",&itype)) itype=0;

    if (!sf_getint("nts",&nts)) sf_error("No nts= in input"); 
    /* maximum time-lag */

    if (!sf_getint("ntw",&ntw)) sf_error("No ntw= in input");
    /* half time-window size */

    /* 2-D section traces contain instantaneous seismic phase */
    trc = sf_floatalloc2(nt,nx);

    sf_floatread(trc[0],nt*nx,in);

    t1 = sf_floatalloc(2*ntw+1);
    t2 = sf_floatalloc(2*ntw+1);

    ntcc = nt - 2*(nts + ntw);

    cc = sf_floatalloc2(ntcc,nx-1);

    for (j = 0; j < 2*ntw+1; j++) {
	t1[j] = 0.0;
	t2[j] = 0.0;
    }

    for (i = 0; i < nx-1; i++) { /* traces loop */

	for (k = nts+ntw; k < nt-nts-ntw; k++) { /* time position loop */

	    t = ot + k*dt;

            /* trace 2 time window */
	    for (j = -ntw; j < ntw+1; j++) t2[j] = trc[i+1][k+j];

	    cm = 0.0; pm = 0.0;

	    for (it = -nts; it < nts+1; it++) { /* time-lags loop */

		ts = it*dt; 

                /* trace 1 time window */
		for (j = 0; j < 2*ntw+1; j++) t1[j] = trc[i][k+j+it];

                /* cross-correlation modulus and phase */
		circ_corr(t1,t2,2*ntw+1,&m,&phi);

                /* extract phase when modulus is maximum */
		if (m > cm) {
		    cm = m; pm = phi;
		}

	    }

            /* coherency is phase of cross-correlation with maximum modulus */
            /* stability with local analysis based on several traces */
	    cc[i][k] = pm;

	}

    }

    sf_floatwrite (cc[0],ntcc*(nx-1),out);
    
    exit (0);
}

/* 	$Id$	 */

