/* Generate angle gathers using the Zoeppritz equation */
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
    int ia, na, icoef, j, i1, i2, n1, n2, n3;
    float vp1,vp2,vs1,vs2,rho1,rho2;
    float a0, da, a, rc[4], ang[4], *r;
    float *vpt, *vst, *rhot;
    sf_file vp, vs, rho, out;
    
    sf_init(argc,argv);
    vp = sf_input("in");
    vs = sf_input("vs");
    rho = sf_input("rho");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(vp)) sf_error("Need float type");
    if (!sf_histint(vp, "n1", &n1)) sf_error("Need n1=");
    n2 = sf_leftsize(vp, 1); /* number of traces */

    vpt = sf_floatalloc(n1);
    vst = sf_floatalloc(n1);
    rhot = sf_floatalloc(n1);
    
    n3 = sf_shiftdim(vp, out, 1);
    if (n3 != n1*n2) sf_error("size mismatch");

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

    
    for (i2=0; i2 < n2; i2++) {
      sf_floatread(vpt,n1,vp);
      sf_floatread(vst,n1,vs);
      sf_floatread(rhot,n1,rho);

      vp1 = vpt[0];
      vs1 = SF_MAX(vst[0],SF_EPS);
      rho1 = rhot[0];
      for (i1=0; i1 < n1; i1++) {
	vp2 = vpt[i1];
	vs2 = SF_MAX(vst[i1],SF_EPS);
	rho2 = rhot[i1];
	
	for (ia=0; ia < na; ia++) {
	  a = a0 + ia*da;
	  a = incp? sinf(a)/vp1: sinf(a)/vs1;

	  zoeppritz (icoef,vp1,vp2,vs1,vs2,rho1,rho2,incp,a,rc,ang);
	  r[ia] = rc[j];
	}

	sf_floatwrite(r,na,out);

	vp1 = vp2;
	vs1 = vs2;
	rho1 = rho2;
      }
    }
      
    exit(0);
}
