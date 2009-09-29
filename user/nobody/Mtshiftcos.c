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

    float v,a,n,f;
    int         fint;

    sf_axa ax,az,av,aa;
    int ix,iz,iv,ia;

    float   **stk, **ang, *vel, *tmp;
    sf_file  Fstk,  Fang, Fvel;

    sf_init (argc,argv);

    /*------------------------------------------------------------*/
    Fstk = sf_input("in");
    Fvel = sf_input("velocity");
    Fang = sf_output("out");

    if (SF_FLOAT != sf_gettype(Fstk)) sf_error("Need float input");

    az=sf_iaxa(Fstk,1); nz=sf_n(az);
    av=sf_iaxa(Fstk,2); nv=sf_n(av);
    ax=sf_iaxa(Fstk,3);

    if (!sf_getint  ("na",&na)) na=nv;       
    if (!sf_getfloat("da",&da)) da=1./(nv-1);
    if (!sf_getfloat("a0",&a0)) a0=0.;         

    aa = sf_maxa(na,a0,da);
    sf_oaxa(Fang,aa,2);

    if (!sf_getint("extend",&ext)) ext=4;       /* tmp extension */
    /*------------------------------------------------------------*/

    stk = sf_floatalloc2(nz,nv);
    ang = sf_floatalloc2(nz,na);
    tmp = sf_floatalloc(nv);
    vel = sf_floatalloc(nz);

    sft = fint1_init(ext,nv,0);
    
    for (ix = 0; ix < ax.n; ix++) {
	sf_floatread(vel   ,nz   ,Fvel);	
	sf_floatread(stk[0],nz*nv,Fstk);
	
	/*------------------------------------------------------------*/
	for (iz = 0; iz < az.n; iz++) {
	    for (iv = 0; iv < av.n; iv++) {
		tmp[iv] = stk[iv][iz];
	    }
	    fint1_set(sft,tmp);
	    v = vel[iz];
	    
	    for (ia=0; ia < aa.n; ia++) {
		a = aa.o+ia*aa.d;      /* ang */
		a = cosf(a/180*SF_PI); /* cos */

		n = v / a;             /* nu = v / cos */
		f = (n - av.o) / av.d;
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
