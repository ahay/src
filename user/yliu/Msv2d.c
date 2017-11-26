/* Velocity and heterogeneity parameter convert to dip. */
/*
  Copyright (C) 2016 Jilin University
  
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
#include "mutter.h"
#include <stdio.h>
#include <math.h>

int main (int argc, char* argv[])
{
    int n1, n2, xn, i, j, k, m, nr;
    float o1, d1, xd, xo, *vel, *dip, t, t1, x2,*ani;
    float t0=0., tp=0., slope0=0., slopep=0., x=0., x1=0., v0=0., *data=NULL;
    bool abs, half, inner, hyper, mute;
    sf_file in, out,anisotropy;
  
    sf_init (argc, argv);

    in = sf_input("in");
    out = sf_output("out");
   
    anisotropy = sf_input("anisotropy");/*anisotropy */
   
    if (!sf_histint(anisotropy,"n1",&nr)) sf_error("No n1= in anisotropy");
  
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("n",&xn)) xn = 32;
    /* offset number */
    if (!sf_getfloat("d",&xd)) xd = 12.5;
    /* offset interval */
    if (!sf_getfloat("o",&xo)) xo = 0.;
    /* offset origin */
    
    if (!sf_getbool("mute",&mute)) mute=false;
    /* if y, use mutter */

    if (!sf_getbool("half",&half)) half=false;
    /* if y, half-offset instead of full offset */

    vel = sf_floatalloc(n1);
    dip = sf_floatalloc(n1*xn);
    
    sf_shiftdim(in, out, 1);
    sf_putint(out,"n2",xn);
    sf_putfloat(out,"d2",xd);
    sf_putfloat(out,"o2",xo);

    if (half) {
	xd *= 2.;
	xo *= 2.;
    }

    if (mute) {
	
	if (!sf_getfloat("tp",&tp)) tp=0.150; 
        /* end time (available when mute=y) */

	if (!sf_getfloat("t0",&t0)) t0=0.; 
	/* starting time (available when mute=y) */

	if (!sf_getfloat("v0",&v0)) v0=10000;  
        /* velocity (available when mute=y) */
	
	if (!sf_getfloat("x0",&x1)) x1=0.;
	/* starting space (available when mute=y) */
	
	if (!sf_getbool("abs",&abs)) abs=true;
	/* if y, use absolute value |x-x0| (available when mute=y) */
	
	if (!sf_getbool("inner",&inner)) inner=false;
	/* if y, do inner muter (available when mute=y) */
	
	if (!sf_getbool("hyper",&hyper)) hyper=false;
	/* if y, do hyperbolic mute (available when mute=y) */
	
	slope0=1./v0;
	slopep=slope0;
	
	if (hyper) {
	    slope0 *= slope0;
	    slopep *= slopep;
	}
	
	data = sf_floatalloc(n1);

        /* Initialize muting parameter */
	mutter_init(n1,o1-t0,d1,abs,inner,hyper);
    }   

    for (k=0; k < n2; k++) {
	sf_warning("slice %d of %d;",k+1,n2);
	sf_floatread(vel,n1,in);

	ani = sf_floatalloc(nr);

	/* Read auxiliary homogeneity parameter S file */
	sf_floatread(ani,nr,anisotropy);

	for (j=0; j < xn; j++) {
	    for (i=0; i < n1; i++)  {
		t = o1+i*d1;
		x2 = xo+j*xd;
		
		m = 0;
		t1 = o1 + m*d1;
		
                /* Find t0 index[variable "m"] for each (t,x) coordinate */
		while(t>(t1*(1-1/ani[m])+sqrt(t1*t1+ani[m]*x2*x2/(vel[m]*vel[m]+FLT_EPSILON))/ani[m])) {
		    m++;
		    t1 = o1 + m*d1;
		}

		/* Calculate sigma(t,x)[variable "dip"] in equ 6 */
		dip[j*n1+i] = x2*xd/
		    (FLT_EPSILON+d1*vel[m]*vel[m]*(ani[m]*(t-t1)+t1));
	    }
	}
	
        /* Muting by using velocity, which is the same as "sfmutter" */
	if (mute) {
	    for (j=0; j < xn; j++) { 
		x = xo + j*xd;	    
		x -= x1;
		if (hyper) x *= x;
		
		for (i=0; i < n1; i++) {
		    data[i] = dip[j*n1+i];
		}
		mutter (tp,slope0,slopep, x, data, false);
		for (i=0; i < n1; i++) {
		    dip[j*n1+i] = data[i];
		}
	    }
	}
	sf_floatwrite(dip,n1*xn,out);
    }

    sf_warning(".");
    exit (0);
}
