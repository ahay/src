/* Compute cos(theta) from 1/|pm| for time-shift imaging condition */
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

#include <math.h>

#include <rsf.h>

#include "fint1.h"

int main (int argc, char* argv[])
{
    fint1 sft;
    int  ext;

    float a,n,f,da,a0;
    int   fint,na,nx,nz,nt;

    sf_axis ax,az,at,aa;
    int ix,iz,it,ia;

    float   **stk, **ang, *tmp;
    sf_file  Fstk,  Fang;

    sf_init (argc,argv);

    /*------------------------------------------------------------*/
    Fstk = sf_input("in");
    Fang = sf_output("out");

    if (SF_FLOAT != sf_gettype(Fstk)) sf_error("Need float input");

    az=sf_iaxa(Fstk,1); nz=sf_n(az);
    at=sf_iaxa(Fstk,2); nt=sf_n(at);
    ax=sf_iaxa(Fstk,3); nx=sf_n(ax);

    if (!sf_getint  ("na",&na)) na=nt;       
    if (!sf_getfloat("da",&da)) da=90/(nt-1);
    if (!sf_getfloat("a0",&a0)) a0=0.;
    aa = sf_maxa(na,a0,da);
    sf_oaxa(Fang,aa,2);

    if (!sf_getint("extend",&ext)) ext=4;       /* tmp extension */
    /*------------------------------------------------------------*/

    stk = sf_floatalloc2(nz,nt);
    ang = sf_floatalloc2(nz,na);
    tmp = sf_floatalloc(nt);

    sft = fint1_init(ext,nt, 0);
    
    for (ix = 0; ix < nx; ix++) {
	sf_floatread(stk[0],nz*nt,Fstk);
	
	/*------------------------------------------------------------*/
	for (iz = 0; iz < nz; iz++) {
	    for (it = 0; it < nt; it++) {
		tmp[it] = stk[it][iz];
	    }
	    fint1_set(sft,tmp);
	    
	    for (ia=0; ia < na; ia++) {
		a = a0+ia*da;      /* ang */
		n = tanf(a/180*SF_PI); /* tan */

		f = (n - sf_o(at)) / sf_d(at);
		fint = f;

		if (fint >= 0 && fint < nt) {
		    ang[ia][iz] = fint1_apply(sft,fint,f-fint,false);
		} else {
		    ang[ia][iz] = 0.;
		}
	    }
	}
	/*------------------------------------------------------------*/
	    
	sf_floatwrite(ang[0],nz*na,Fang);
    }
	
    exit (0);
}
