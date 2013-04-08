/* Analytic escape solutions in phase space for constant gradient of slowness squared */
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
#include <assert.h>

static int nz,nx,na;
static float oz,ox,oa;
static float dz,dx,da;

void solve_exit(float sexit,         /* exit position to solve */
                float s,             /* initial position */
                float ps,            /* initial slowness component */
                float gs,            /* gradient component */
                float *sol1,         /* root 1 of quadratic equation */
                float *sol2          /* root 2 of quadratic equation */
                ) 
{
    float sm1,sm2,b,c,d,q;
    
    gs *= 0.5;
    sm1 = sm2 = -1e6f;

    if (fabs(gs) < 1e-6f) {

	if (fabs(ps) > 1e-6f) {
	    sm1 = (sexit-s)/ps;
	} else { 
	    assert(0);   
	}  

    } else {

	/* sm^2 + b*sm + c = 0 */

	b = ps/gs;
	c = (s-sexit)/gs;

	if (fabs(b) < 1e-6f) {

	    d = -c;
	    if (d > 1e-6f) sm1 = sqrt(d);

	} else {

	    d = b*b - 4.0 * c;

	    if (fabs(d) < 1e-6f) {

		sm1 = -0.5*b;

	    } else if (d > 1e-6f) {

		q = -0.5*(b + SF_SIG(b)*sqrt(d));

		sm1 = q;
		sm2 = c/q;
	    
	    }
	}
    }

    *sol1 = sm1; 
    *sol2 = sm2;
    return;
}


int main(int argc, char* argv[])
{
    int iq;                        /* escape variables switch */
    int iz,ix,ia,is;

    float gx,gz,sc,xc,zc;          /* slowness parameters */
    float xmax,zmax;
    float a,x,z,s,px,pz;

    float ss[8],sm,smp;
    float t_min,t_tmp;

    float *t;                    /* escape variable */

    sf_file in,out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("iq",&iq)) sf_error("Need iq=");
    /* switch for escape variable 0=x, 1=a, 2=t, 3=z */

    if (!sf_getfloat("gx",&gx)) sf_error("Need gx=");
    /* x-gradient of slowness^2 */

    if (!sf_getfloat("gz",&gz)) sf_error("Need gz=");
    /* z-gradient of slowness^2 */

    if (!sf_getfloat("sc",&sc)) sf_error("Need sc=");
    /* slowness constant */

    if (!sf_getfloat("xc",&xc)) sf_error("Need xc=");
    /* x reference */

    if (!sf_getfloat("zc",&zc)) sf_error("Need zc=");
    /* z reference */

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");

    if (!sf_histint(in,"n3",&na)) sf_error("No n1=");
    if (!sf_histfloat(in,"d3",&da)) sf_error("No d1=");
    if (!sf_histfloat(in,"o3",&oa)) sf_error("No o1=");

    if (!sf_histint(in,"n2",&nx)) sf_error("No n2=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2=");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n3=");
    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d3=");
    if (!sf_histfloat(in,"o1",&oz)) sf_error("No o3=");
    /* angle in degrees from -180 to 180 */

    /* memory allocations */
    t = sf_floatalloc(nz);

    /* convert to radians */
    oa *= SF_PI/180.;
    da *= SF_PI/180.;

    /* Boundaries */
    xmax = ox + (nx-1)*dx;
    zmax = oz + (nz-1)*dz;

    for (ia = 0; ia < na; ia++) {
			
	a = oa + ia*da;

	for (ix = 0; ix < nx; ix ++) {

	    x = ox + ix*dx;

	    for (iz = 0; iz < nz; iz++) {
 
		z = oz + iz*dz;
		s = sqrt(sc*sc + 2.*(gx*(x-xc) + gz*(z-zc)));

		px = -s*sin(a);
		pz = -s*cos(a);

		if ((iz == 0    && pz < -1e-6) ||  
                    (iz == nz-1 && pz >  1e-6) ||  
                    (ix == 0    && px < -1e-6) || 
                    (ix == nx-1 && px >  1e-6)) {
                    /* x and z boundaries */

		    switch(iq) {
			case 0:
                            /* escape location x */
			    t[iz] = x;
			    break;
			case 1:
                            /* escape angle a */
			    t[iz] = a*180./SF_PI;
			    break;
			case 2:
                            /* escape traveltime */
			    t[iz] = 0.0;
			    break;
			case 3:
                            /* escape depth z */
			    t[iz] = z;
			    break;
		    }

		} else {

		    solve_exit(ox,x,px,gx,ss,ss+1);
		    solve_exit(xmax,x,px,gx,ss+2,ss+3);

		    solve_exit(oz,z,pz,gz,ss+4,ss+5);
		    solve_exit(zmax,z,pz,gz,ss+6,ss+7);

		    /* select lowest positive traveltime */

		    sm = 1e6f;
		    t_min = 1e6f;

		    for (is = 0; is < 8; is++) {

			smp = ss[is];
			t_tmp =  (s*s*smp + (px*gx + pz*gz)*smp*smp + (gx*gx + gz*gz)*smp*smp*smp/3.f);
			
			if (t_tmp > 1e-6 && t_tmp < t_min ) {
			    sm = smp;
			    t_min = t_tmp;
			}

		    }

		    switch(iq) {
			case 0:
                            /* escape location x */
			    t[iz] = x + px*sm + 0.5*gx*sm*sm;
			    break;
			case 1:
                            /* escape angle a */
			    t[iz] = atan2(-(px+gx*sm),-(pz+gz*sm))*180./SF_PI;
			    break;
			case 2:
                            /* escape traveltime */
			    t[iz] = s*s*sm + (px*gx + pz*gz)*sm*sm + 0.3333333333*(gx*gx + gz*gz)*sm*sm*sm;
			    break;
			case 3:
                            /* escape depth z */
			    t[iz] = z + pz*sm + 0.5*gz*sm*sm;
			    break;
		    }

		} /* if */ 
		
	    } /* iz */

	    
	    /* output */
	    sf_floatwrite(t,nz,out);

	} /* ix */

    } /* ia */
    
    exit(0);
}

