/* 2-D synthetic data of conflicting dips. */
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

#include "lint1.h"
#include "random.h"

#define NT1 15
#define NT2 25

int main (int argc, char* argv[])
{
    int nt, nx, it, ix;
    float **data, t[NT2], amp[NT2], t0, x0, dt, dx;
    sf_file mod;

    sf_init (argc,argv);
    mod = sf_output("out");
    sf_setformat(mod,"native_float");

    if (!sf_getint("n1",&nt)) nt=150;
    if (!sf_getint("n2",&nx)) nx=80;
    sf_putint(mod,"n1",nt);
    sf_putint(mod,"n2",nx);    

    t0 = 0.; 
    x0= 0.;
    dt = 1./nt; 
    dx= 1.;

    sf_putfloat(mod,"d1",dt);
    sf_putfloat(mod,"o1",t0);
    sf_putfloat(mod,"d2",dx);
    sf_putfloat(mod,"o2",x0);

    data = sf_floatalloc2(nt,nx);

    random_init (1992);
    
    for (ix=0; ix < nx; ix++) {
	for (it=0; it < nt; it++) {
	    data[ix][it] = 1.e-2*(random0()-0.5);
	}
    }

    for (ix=0; ix < nx; ix++) {
	for (it=0; it < NT1; it++) {
	    amp[it] = (ix+1.)/ nx;
	    t[it] = 0.4  * (ix+1.) / nx + (-2*nt + it*(nt/3))*dt;
	}
	lint1_init (t0,dt,t);
	lint1_lop (true,true,nt,NT1,data[ix],amp);
    }

    for (ix=0; ix < nx; ix++) {
	for (it=0; it < NT2; it++) {
	    amp[it] = (float) (nx-ix)/ nx;
	    t[it]   = 0.7  * (nx-ix-1.) / nx + (-2*nt + it*(nt/5))*dt;
	}
	lint1_init (t0,dt,t);
	lint1_lop (true,true,nt,NT2,data[ix],amp);
    }

    sf_floatwrite(data[0],nt*nx,mod);

    exit(0);
}

/* 	$Id$	 */
