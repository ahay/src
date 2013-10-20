/* Convert Dix velocity to interval velocity (2-D). */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

int main(int argc, char *argv[])
{
    bool verb;
    int ix, it, nx, nt, niter, rect;
    float ox, dx, dt, norm;
    float **v, **vi, *x0, *z0, *x1, *z1, *x2, *num, *den;
    sf_file vdix, zx;

    sf_init(argc,argv);
    vdix = sf_input("in");
    zx = sf_output("out");

    if (!sf_histint(vdix,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(vdix,"n2",&nt)) sf_error("No n2= in input");

    if (!sf_histfloat(vdix,"o1",&ox)) sf_error("No o2= in input");
    if (!sf_histfloat(vdix,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(vdix,"d2",&dt)) sf_error("No d2= in input");

    sf_shiftdim(vdix,zx,2);
    sf_putint(zx,"n2",2);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("rect",&rect)) rect=3;
    /* lateral smoothing */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    v  = sf_floatalloc2(nx,nt);
    vi = sf_floatalloc2(nx,nt);

    x0 = sf_floatalloc(nx);
    z0 = sf_floatalloc(nx);
    x1 = sf_floatalloc(nx);
    z1 = sf_floatalloc(nx);
    x2 = sf_floatalloc(nx);

    sf_floatread(v[0],nx*nt,vdix);
    
    for (ix=0; ix < nx; ix++) {
	x0[ix] = 0.;
	z0[ix] = 0.;
	
	x1[ix] = 0.;
	z1[ix] = dt*v[0][ix];

	x2[ix] = ox+ix*dx;
    }

    sf_floatwrite(z0,nx,zx);
    sf_floatwrite(x2,nx,zx);

    sf_floatwrite(z1,nx,zx);
    sf_floatwrite(x2,nx,zx);

    /* dimensionless velocity */
    for (it=0; it < nt; it++) {
	for (ix=0; ix < nx; ix++) {
	    v[it][ix] *= dt/dx;
	    vi[it][ix] = 1.0/v[it][ix];
	}
    }

    sf_divn_init(1,nx,&nx,&rect,niter,verb);

    den = sf_floatalloc(nx);
    num = sf_floatalloc(nx);

    for (it=2; it < nt; it++) {
	norm = 0.;
	for (ix=0; ix < nx; ix++) {
	    den[ix] = vi[it-1][ix]+0.25*(vi[it][ix]-vi[it-2][ix]);
	    norm += den[ix]*den[ix];
	}
	norm = sqrtf(nx/norm);

	for (ix=0; ix < nx; ix++) {
	    den[ix] *= norm;
	    if (0==ix) {
		num[ix] = -norm*(v[it-1][ix]*(z1[ix]-2*z1[ix+1]+z1[ix+2])+
				 vi[it-1][ix]*(z0[ix]-2*z1[ix])+
				 (v[it-1][ix+1]-v[it-1][ix])*(z1[ix+1]-z1[ix])-
				 0.25*(vi[it][ix]-vi[it-2][ix])*z0[ix]);
	    } else if (nx-1==ix) {
		num[ix] = -norm*(v[it-1][ix]*(z1[ix-2]-2*z1[ix-1]+z1[ix])+
				 vi[it-1][ix]*(z0[ix]-2*z1[ix])+
				 (v[it-1][ix]-v[it-1][ix-1])*(z1[ix]-z1[ix-1])-
				 0.25*(vi[it][ix]-vi[it-2][ix])*z0[ix]);
	    } else {
		num[ix] = -norm*(v[it-1][ix]*(z1[ix-1]-2*z1[ix]+z1[ix+1])+
				 vi[it-1][ix]*(z0[ix]-2*z1[ix])+
				 0.25*(v[it-1][ix+1]-v[it-1][ix-1])*(z1[ix+1]-z1[ix-1])-
				 0.25*(vi[it][ix]-vi[it-2][ix])*z0[ix]);
	    }
	    z0[ix] = z1[ix];
	}
	sf_divn(num,den,z1);
	sf_floatwrite(z1,nx,zx);  

	for (ix=0; ix < nx; ix++) {
	    if (0==ix) {
		num[ix] = -norm*(v[it-1][ix]*(x1[ix]-2*x1[ix+1]+x1[ix+2])+
				 vi[it-1][ix]*(x0[ix]-2*x1[ix])+
				 (v[it-1][ix+1]-v[it-1][ix])*(x1[ix+1]-x1[ix])-
				 0.25*(vi[it][ix]-vi[it-2][ix])*x0[ix]+
				 (v[it-1][ix+1]-v[it-1][ix])*dx);
	    } else if (nx-1==ix) {
		num[ix] = -norm*(v[it-1][ix]*(x1[ix-2]-2*x1[ix-1]+x1[ix])+
				 vi[it-1][ix]*(x0[ix]-2*x1[ix])+
				 (v[it-1][ix]-v[it-1][ix-1])*(x1[ix]-x1[ix-1])-
				 0.25*(vi[it][ix]-vi[it-2][ix])*x0[ix]+
				 (v[it-1][ix]-v[it-1][ix-1])*dx);
	    } else {
		num[ix] = -norm*(v[it-1][ix]*(x1[ix-1]-2*x1[ix]+x1[ix+1])+
				 vi[it-1][ix]*(x0[ix]-2*x1[ix])+
				 0.25*(v[it-1][ix+1]-v[it-1][ix-1])*(x1[ix+1]-x1[ix-1])-
				 0.25*(vi[it][ix]-vi[it-2][ix])*x0[ix]+
				 0.5*(v[it-1][ix+1]-v[it-1][ix-1])*dx);
	    }
	    x0[ix] = x1[ix];
	}
	sf_divn(num,den,x1);

	for (ix=0; ix < nx; ix++) {
	    x2[ix] = x1[ix] + ox+ix*dx;
	}
	sf_floatwrite(x2,nx,zx);
    }
    
    exit(0);
}
