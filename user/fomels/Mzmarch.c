/* Phase-space traveltime (marching in z) */
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

static float deltat(int scheme, float dz2, 
		    float dz, float dx, 
		    float s, float s2, 
		    float a, float a2, 
		    const float* grad, const float* grad2);

int main(int argc, char* argv[])
{
    bool verb;
    char *what;
    int nz, nx, na, iz, ix, ia, order, jx, ja, scheme;
    float x, dx, x0, a, da, a0, dz, dz2, t, sx, sa;
    float amax, xmax, x2, a1, a2, s2, tn;
    float **next, **vxz, s, grad[2], grad2[2];
    sf_eno2 prev, slow;
    sf_file vel, out;

    sf_init(argc,argv);
    vel = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2=");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(vel,"o2",&x0)) sf_error("No o2=");

    if (!sf_histint(vel,"n1",&nz)) sf_error("No n2=");
    if (!sf_histfloat(vel,"d1",&dz)) sf_error("No d2=");

    if (!sf_getint("na",&na)) sf_error("Need na=");
    /* angle samples */
    if (!sf_getfloat("da",&da)) sf_error("Need da=");
    /* angle sampling */
    if (!sf_getfloat("a0",&a0)) sf_error("Need a0=");
    /* starting velocity */

    sf_putint(out,"n1",na);
    sf_putint(out,"n2",nx);
    sf_putint(out,"n3",nz);
    sf_putfloat(out,"d1",da);
    sf_putfloat(out,"d2",dx);
    sf_putfloat(out,"d3",dz);
    sf_putfloat(out,"o1",a0);
    sf_putfloat(out,"o2",x0);
    sf_putfloat(out,"o3",0.);
    sf_putstring(out,"label1","Angle");
    sf_putstring(out,"unit1","\\^o\\_");	

    /* convert to radians */
    a0 *= SF_PI/180.;
    da *= SF_PI/180.;

    amax = a0 + (na-1)*da;
    xmax = x0 + (nx-1)*dx;

    if (!sf_getint("order",&order)) order=3;
    /* interpolation order */

    if (!sf_getint("scheme",&scheme)) scheme=2;
    /* finite-difference order */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */

    what = sf_getstring("what");
    /* what to compute (t,x,z,a) */

    if (NULL == what) {
	what = "time";
    } else {
	if (what[0] != 't' &&
	    what[0] != 'x' &&
	    what[0] != 'z' &&
	    what[0] != 'a') 
	    sf_error("Need what=t|x|z|a");
    }

    slow = sf_eno2_init (order,nz,nx);
    vxz = sf_floatalloc2(nz,nx);

    sf_floatread(vxz[0],nz*nx,vel);
    /* convert velocity to slowness */
    for (ix=0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
	    vxz[ix][iz] = 1./vxz[ix][iz];
	}
    }
    sf_eno2_set(slow,vxz);
    free(*vxz);
    free(vxz);

    prev = sf_eno2_init(order,na,nx);
    next = sf_floatalloc2(na,nx);
    for (ix=0; ix < nx; ix++) {
	x = x0+ix*dx;
	for (ia=0; ia < na; ia++) {
	    switch (what[0]) {
		case 't':
		case 'z':
		    next[ix][ia] = 0.;
		    break;
		case 'x':
		    next[ix][ia] = x;
		    break;
		case 'a':
		    a = a0+ia*da;
		    next[ix][ia] = a*180./SF_PI;
	    }
	}
    }
    sf_eno2_set(prev,next);
    sf_floatwrite(next[0],na*nx,out);

    for (iz=1; iz < nz; iz++) {
	if (verb) sf_warning("depth %d of %d",iz,nz);
	for (ix=0; ix < nx; ix++) {	
	    x = x0+ix*dx;

	    /* slowness and gradient at starting point */
	    sf_eno2_apply (slow,iz,ix,0.,0.,&s,grad,BOTH);

	    for (ia=0; ia < na; ia++) {
		a = a0+ia*da;

		if (scheme > 1) {
		    /* half a step in a */
		    a1 = a + 0.5*(grad[1]/dx+
				  grad[0]*tanf(a)/dz)*dz/s;
	
		    /* check if angle is in range */
		    if (a1 < a0) {
			a1=a0;
		    } else if (a1 > amax) {
			a1=amax;
		    } 
		} else {
		    a1 = a;
		}

		tn = tanf(a1);

		/* full step in x */
		x2 = x + tn*dz;
		
		/* check if space is in range */
		if (x2 < x0) {
		    x2 = x0;
		    jx = 0;
		    sx = 0.;
		    dz2 = (x0-x)/tn;
		} else if (x2 > xmax) {
		    x2 = xmax;
		    jx = nx-1;
		    sx = 0.;
		    dz2 = (xmax-x)/tn;
		} else {
		    sx = (x2-x0)/dx; 
		    jx=floorf(sx);
		    sx -= jx;
		    dz2 = dz;
		}

		/* slowness and gradient at new location */
		sf_eno2_apply (slow,iz-1,jx,1.-dz2/dz,sx,
			       &s2,grad2,BOTH);
		
		if (scheme > 1) {
		    /* half a step in a */
		    a2 = a1 + 0.5*(grad2[1]/dx+
				   grad2[0]*tn/dz)*dz2/s2;
		} else {
		    /* full step in a */
		    a2 = a1 + (grad2[1]/dx+
			       grad2[0]*tn/dz)*dz2/s2;
		}

		/* check if angle is in range */
		if (a2 < a0) {
		    a2=a0;
		} else if (a2 > amax) {
		    a2=amax;
		} 

		if (dz2 != dz) { /* exit from the side */
		    switch (what[0]) {
			case 't':			    
			    t = deltat(scheme,dz,dz,dx,
				       s,s2,a,a2,grad,grad2);
			    break;
			case 'x':
			    t = x2;
			    break;
			case 'z':
			    t = iz*dz-dz2;
			    break;
			case 'a':
			    t = a2*180./SF_PI;
			    break;
		    }
		} else { /* exit to the previous level */
		    sa = (a2-a0)/da; 
		    ja=floorf(sa);
		    sa -= ja;

		    sf_eno2_apply (prev,ja,jx,sa,sx,
				   &t,NULL,FUNC);
		    if ('t'==what[0]) 
			t += deltat(scheme,dz,dz,dx,
				    s,s2,a,a2,grad,grad2);
		}
		
		next[ix][ia] = t;
	    }
	}
	sf_eno2_set(prev,next);
	sf_floatwrite(next[0],na*nx,out);
    }
    
    exit(0);
}

static float deltat(int scheme, float dz2, 
		    float dz, float dx, 
		    float s, float s2, 
		    float a, float a2, 
		    const float* grad, const float* grad2)
{
    float t, t2, c, c2;

    c  = 1/cosf(a );
    c2 = 1/cosf(a2);

    t = 0.5*(s*c+s2*c2);
    if (scheme < 2) return t*dz2;

    c  *= c *c ;
    c2 *= c2*c2;

    t2 = t + ((grad [1]*sinf(2*a )/dx+
	       grad [0]*cosf(2*a )/dz)*c-
	      (grad2[1]*sinf(2*a2)/dx+
	       grad2[0]*cosf(2*a2)/dz)*c2)*dz2/12.0;
    if (t2 > 0.) return t2*dz2;
    return t;
}
