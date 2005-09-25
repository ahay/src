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

    float a,n,f;
    int       fint;

    axa ax,az,at,aa;
    int ix,iz,it,ia;

    float   **stk, **ang, *tmp;
    sf_file  Fstk,  Fang;

    sf_init (argc,argv);

    /*------------------------------------------------------------*/
    Fstk = sf_input("in");
    Fang = sf_output("out");

    if (SF_FLOAT != sf_gettype(Fstk)) sf_error("Need float input");

    iaxa(Fstk,&az,1);
    iaxa(Fstk,&at,2);
    iaxa(Fstk,&ax,3);

    if (!sf_getint  ("na",&aa.n)) aa.n=    at.n;       
    if (!sf_getfloat("da",&aa.d)) aa.d=90/(at.n-1);
    if (!sf_getfloat("a0",&aa.o)) aa.o=0.;         
    oaxa(Fang,&aa,2);

    if (!sf_getint("extend",&ext)) ext=4;       /* tmp extension */
    /*------------------------------------------------------------*/

    stk = sf_floatalloc2(az.n,at.n);
    ang = sf_floatalloc2(az.n,aa.n);
    tmp = sf_floatalloc(      at.n);

    sft = fint1_init(ext, at.n, 0);
    
    for (ix = 0; ix < ax.n; ix++) {
	sf_floatread(stk[0],az.n*at.n,Fstk);
	
	/*------------------------------------------------------------*/
	for (iz = 0; iz < az.n; iz++) {
	    for (it = 0; it < at.n; it++) {
		tmp[it] = stk[it][iz];
	    }
	    fint1_set(sft,tmp);
	    
	    for (ia=0; ia < aa.n; ia++) {
		a = aa.o+ia*aa.d;      /* ang */
		n = tanf(a/180*SF_PI); /* tan */

		f = (n - at.o) / at.d;
		fint = f;

		if (fint >= 0 && fint < at.n) {
		    ang[ia][iz] = fint1_apply(sft,fint,f-fint,false);
		} else {
		    ang[ia][iz] = 0.;
		}
	    }
	}
	/*------------------------------------------------------------*/
	    
	sf_floatwrite(ang[0],az.n*aa.n,Fang);
    }
	
    exit (0);
}
