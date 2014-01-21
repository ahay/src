/* Testing Zoeppritz equation */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "zoeppritz.h"

int main(int argc, char* argv[])
{
    bool incp, outp, refl;
    int ia, na, icoef, j;
    float vp1,vp2,vs1,vs2,rho1,rho2;
    float a0, da, a, rc[4], ang[4], *r;
    sf_file out;
    
    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");

    if (!sf_getint("na",&na)) na=90;
    /* number of angles */

    if (!sf_getfloat("a0",&a0)) a0=0.;
    /* first angle */
 
    if (!sf_getfloat("da",&da)) da=90./na;
    /* angle increment */
    
    sf_putint(out,"n1",na);
    sf_putfloat(out,"o1",a0);
    sf_putfloat(out,"d1",da);

    sf_putstring(out,"label1","Incident Angle");
    sf_putstring(out,"unit1","\\^o\\_");

    a0 *= SF_PI/180.;
    da *= SF_PI/180.;

    r = sf_floatalloc(na);

    if (!sf_getint("icoef",&icoef)) icoef=4;
    /* [1,2,3,4] particle displacement, displacement potential, energy, real part */

    if (!sf_getfloat("vp1",&vp1)) sf_error("Need vp1=");
    if (!sf_getfloat("vp2",&vp2)) sf_error("Need vp2=");
    if (!sf_getfloat("vs1",&vs1)) sf_error("Need vs1=");
    if (!sf_getfloat("vs2",&vs2)) sf_error("Need vs2=");
    if (!sf_getfloat("rho1",&rho1)) rho1=1.;
    if (!sf_getfloat("rho2",&rho2)) rho2=1.;
    
    if (!sf_getbool("incp",&incp)) incp=true;
    /* incident P (or S) */

    if (!sf_getbool("outp",&outp)) outp=true;
    /* rellected/transmitted P (or S) */

    if (!sf_getbool("refl",&refl)) refl=true;
    /* reflection or transmission */


    if (outp) {
	j = refl? 0:2;
    } else {
	j = refl? 1:3;
    }

    for (ia=0; ia < na; ia++) {
	a = a0 + ia*da;
	a = incp? sinf(a)/vp1: sinf(a)/vs1;

	zoeppritz (icoef,vp1,vp2,vs1,vs2,rho1,rho2,incp,a,rc,ang);
	r[ia] = rc[j];
    }

    sf_floatwrite(r,na,out);
    exit(0);
}
