/* Riemannian Wavefield Extrapolation: modeling */
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
#include "rweone.h"

int main(int argc, char* argv[])
{

    sf_file Fi, Fo, Fm, Fr;
    axa ax,az,aw,ar;
    int ix,iz,   ir;

    int method;
    bool verb;

    complex float **dat;
    float         **img;
    float         **aa,**bb,**mm;

    complex float **ab;
    float         **a0,**b0;

    sf_init(argc,argv);

    Fi = sf_input("in");
    Fm = sf_input("abm");
    Fr = sf_input("abr");
    Fo = sf_output("out");

    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getint("method",&method)) method=0;

    iaxa(Fi,&ax,1);       /* x='position' (can be angle) */
    iaxa(Fi,&aw,2);       /* w=frequency */
    iaxa(Fm,&az,1);       /* z=extrapolation (can be time) */
    iaxa(Fr,&ar,2);       /* a,b reference */
    if(method==0) ar.n=1; /* pure F-D */

    if(verb) {
	raxa(ax);
	raxa(az);
	raxa(aw);
	raxa(ar);
    }

    oaxa(Fo,&az,1);
    oaxa(Fo,&ax,2);
    sf_settype(Fo,SF_FLOAT);

    /* read data */
    sf_warning("read data");
    dat = sf_complexalloc2(ax.n,aw.n);
    sf_complexread(dat[0],ax.n*aw.n,Fi);

    /* read ABM */
    sf_warning("read ABM");
    aa = sf_floatalloc2  (az.n,ax.n);
    bb = sf_floatalloc2  (az.n,ax.n);
    mm = sf_floatalloc2  (az.n,ax.n);

    sf_floatread(aa[0],az.n*ax.n,Fm); /* a coef */
    sf_floatread(bb[0],az.n*ax.n,Fm); /* b coef */
    sf_floatread(mm[0],az.n*ax.n,Fm); /* mask */
    sf_warning("read ABM ok");

    /* read ABr */
    sf_warning("read ABr");
    ab = sf_complexalloc2(az.n,ar.n);
    a0 = sf_floatalloc2  (az.n,ar.n);
    b0 = sf_floatalloc2  (az.n,ar.n);

    sf_complexread(ab[0],az.n*ar.n,Fr);
    for(ir=0;ir<ar.n;ir++) {
	for(iz=0;iz<az.n;iz++) {
	    a0[ir][iz] = crealf(ab[ir][iz]);
	    b0[ir][iz] = cimagf(ab[ir][iz]);
	}
    }
    sf_warning("read ABr ok");

    /* allocate image */
    sf_warning("allocate image");
    img = sf_floatalloc2  (az.n,ax.n);
    for(ix=0;ix<ax.n;ix++) {
	for(iz=0;iz<az.n;iz++) {
	    img[ix][iz] = 0.;
	}
    }

    /* model */
    sf_warning("init modeling");
    rweone_init(ax,az,aw,ar,method);

    sf_warning("run modeling");
    rweone_main(dat,img,aa,bb,mm,a0,b0);

    /* write image */
    sf_warning("write image");
    sf_floatwrite(img[0],ax.n*az.n,Fo);

    sf_warning("OK");
}
