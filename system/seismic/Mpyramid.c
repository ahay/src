/* Pyramid transform */
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

int main(int argc, char* argv[])
{
    sf_map4 mo;
    bool inv;
    int i, nx, nu, iw, nw;
    float u0, du, x0, dx, eps;
    sf_complex *ctrace, *ctrace2;
    float *rtrace, *itrace, *str, *rtrace2, *itrace2;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* inversion flag */

    if (inv) {
	if (!sf_histint(in,"n1",&nu)) sf_error("No nu= in input");
	if (!sf_histfloat(in,"d1",&du)) du=1.;
	if (!sf_histfloat(in,"o1",&u0)) u0=0.;

	if (!sf_getint("nx",&nx) && !sf_histint(in,"pyr_nx",&nx)) sf_error("Need nx=");
	if (!sf_getfloat("dx",&dx) && !sf_histfloat(in,"pyr_dx",&dx)) dx=du;
	if (!sf_getfloat("x0",&x0) && !sf_histfloat(in,"pyr_x0",&x0)) x0=u0;
	
	sf_putint(out,"n1",nx);
	sf_putfloat(out,"d1",dx);
	sf_putfloat(out,"o1",x0);
    } else {
	if (!sf_histint(in,"n1",&nx)) sf_error("No nx= in input");
	if (!sf_histfloat(in,"d1",&dx)) dx=1.;
	if (!sf_histfloat(in,"o1",&x0)) x0=0.;

	if (!sf_getint("nu",&nu)) sf_error("Need nu=");
	if (!sf_getfloat("du",&du)) du=dx;
	if (!sf_getfloat("u0",&u0)) u0=x0;
	
	sf_putint(out,"n1",nu);
	sf_putfloat(out,"d1",du);
	sf_putfloat(out,"o1",u0);

	sf_putint(out,"pyr_nx",nx);
	sf_putfloat(out,"pyr_dx",dx);
	sf_putfloat(out,"pyr_x0",x0);
    }

    nw = sf_leftsize(in,1);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    rtrace = sf_floatalloc(nx);
    itrace = sf_floatalloc(nx);
    ctrace = sf_complexalloc(nx);

    str = sf_floatalloc(nx);

    rtrace2 = sf_floatalloc(nu);
    itrace2 = sf_floatalloc(nu);
    ctrace2 = sf_complexalloc(nu);

    mo = sf_stretch4_init (nu, u0, du, nx, eps);

    for (iw=0; iw < nw; iw++) {
	for (i=0; i < nx; i++) {
	    str[i] = iw*(x0+i*dx);
	}
	sf_stretch4_define (mo,str,false);

	if (inv) {
	    sf_complexread(ctrace2,nu,in);
	    
	    for (i=0; i < nu; i++) {
		rtrace2[i] = crealf(ctrace2[i]);
		itrace2[i] = cimagf(ctrace2[i]);
	    }

	    sf_stretch4_invert (false,mo,rtrace,rtrace2);
	    sf_stretch4_invert (false,mo,itrace,itrace2);

	    for (i=0; i < nx; i++) {
		ctrace[i] = sf_cmplx(rtrace[i],itrace[i]);
	    }

	    sf_complexwrite(ctrace,nx,out);
	} else {
	    sf_complexread(ctrace,nx,in);

	    for (i=0; i < nx; i++) {
		rtrace[i] = crealf(ctrace[i]);
		itrace[i] = cimagf(ctrace[i]);
	    }

	    sf_stretch4_apply (false,mo,rtrace,rtrace2);
	    sf_stretch4_apply (false,mo,itrace,itrace2);

	    for (i=0; i < nu; i++) {
		ctrace2[i] = sf_cmplx(rtrace2[i],itrace2[i]);
	    }

	    sf_complexwrite(ctrace2,nu,out);
	} 
    }
    
    exit(0);
}
