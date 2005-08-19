/* Compute shift from pseudo-v to pseudo-tan(theta) */
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

    float v,a,n,f;
    int         fint;

    axa ax,az,av,aa;
    int ix,iz,iv,ia;

    float   **stk, **ang, *vel, *tmp;
    sf_file  Fstk,  Fang, Fvel;

    sf_init (argc,argv);

    /*------------------------------------------------------------*/
    Fstk = sf_input("in");
    Fvel = sf_input("velocity");
    Fang = sf_output("out");

    if (SF_FLOAT != sf_gettype(Fstk)) sf_error("Need float input");

    iaxa(Fstk,&az,1);
    iaxa(Fstk,&av,2);
    iaxa(Fstk,&ax,3);

    if (!sf_getint  ("na",&aa.n)) aa.n=    av.n;       
    if (!sf_getfloat("da",&aa.d)) aa.d=1./(av.n-1);
    if (!sf_getfloat("a0",&aa.o)) aa.o=0.;         

    oaxa(Fang,&aa,2);

    if (!sf_getint("extend",&ext)) ext=4;       /* tmp extension */
    /*------------------------------------------------------------*/

    stk = sf_floatalloc2(az.n,av.n);
    ang = sf_floatalloc2(az.n,aa.n);
    tmp = sf_floatalloc(      av.n);
    vel = sf_floatalloc(az.n      );

    sft = fint1_init(ext, av.n);
    
    for (ix = 0; ix < ax.n; ix++) {
	sf_floatread(vel   ,az.n     ,Fvel);	
	sf_floatread(stk[0],az.n*av.n,Fstk);
	
	/*------------------------------------------------------------*/
	for (iz = 0; iz < az.n; iz++) {
	    for (iv = 0; iv < av.n; iv++) {
		tmp[iv] = stk[iv][iz];
	    }
	    fint1_set(sft,tmp);
	    v = vel[iz];
	    
	    for (ia=0; ia < aa.n; ia++) {
		a = aa.o+ia*aa.d;      /*                    tan     */
		n = v * hypotf(a,1.);  /* nu = v * sqrt( 1 - tan^2 ) */

/*		f = (n - av.o) / av.d;*/
/*		if( a>0. ) {*/
/*		    f = ( v-av.o + SF_ABS(n-v) ) / av.d;*/
/*		} else {*/
/*		    f = ( v-av.o - SF_ABS(n-v) ) / av.d;*/
/*		}*/
		f =  ( v-av.o + SF_SIG(a) * SF_ABS(n-v) ) / av.d;

		fint = f;

		if (fint >= 0 && fint < av.n) {
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
