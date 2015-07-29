/* Generate Lorenz attractor. */
/*
  Copyright (C) 2011 Jilin University
   
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

int main(int argc, char* argv[]) 
{
    int i, niter, n;
    double rho, sigma, beta, x0, y0, z0, x, y, z, dt;
    sf_complex *set;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_complex");

    if (!sf_getint("niter",&niter)) niter=1000;
    /* number of iterations */
    if (!sf_getint("n",&n)) n=niter;
    /* set maximum */
    if (niter > n || niter <= 0 || n <= 0) sf_error("need 0 < niter <= n");
    if (!sf_getdouble("rho",&rho)) rho=28.00;
    /* Rayleigh number */
    if (!sf_getdouble("sigma",&sigma)) sigma=10.00;
    /* Prandtl number */
    if (!sf_getdouble("beta",&beta)) beta=8.00/3.00;
    /* Beta reference */

    if (!sf_getdouble("x0",&x0)) x0=3.051522;
    /* intial x coordinate */
    if (!sf_getdouble("y0",&y0)) y0=1.582542;
    /* intial y coordinate */
    if (!sf_getdouble("z0",&z0)) z0=15.62388;
    /* intial z coordinate */

    if (!sf_getdouble("dt",&dt)) dt=0.0001;
    /* time step */

    sf_putint(out,"n1",n);
    sf_putint(out,"n2",3);

    sf_putfloat(out,"o1",0.);
    sf_putfloat(out,"o2",0.);

    sf_putfloat(out,"d1",1.);
    sf_putfloat(out,"d2",1.);

    set = sf_complexalloc(3*n);

    x = x0; y = y0; z = z0;
    for (i=0; i < niter; i++) {
	set[i] = sf_cmplx(x,y);
	set[n+i] = sf_cmplx(x,z);
	set[2*n+i] = sf_cmplx(y,z);
	x += sigma*(y-x)*dt;
	y += (rho*x-x*z-y)*dt;
	z += (x*y-beta*z)*dt;
    }

    for (i=niter; i < n; i++) {
	set[i] = set[i-1];
	set[n+i] = set[n+i-1];
	set[2*n+i] = set[2*n+i-1];
    }
    sf_complexwrite(set,3*n,out);
 
    exit(0);
}

/* 	$Id: Mlorenz.c 9567 2012-10-29 20:38:18Z sfomel $	 */
