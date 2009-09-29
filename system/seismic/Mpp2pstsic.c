/* Compute angle gathers for time-shift imaging condition */
/*
  Copyright (C) 2006 Colorado School of Mines

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
    int  ext, nx, nz, nv, na;

    float v,a,c,d,n,e,f,g, a0,da, v0,dv;
    int         fint;

    sf_axis ax,az,av,aa;
    int ix,iz,iv,ia;

    float   **stk=NULL, **ang=NULL, *vel=NULL, *dzdx=NULL, *gam=NULL, *tmp=NULL;
    sf_file  Fstk=NULL,  Fang=NULL, Fvel=NULL, Fdip=NULL,  Fgam=NULL;

    sf_init (argc,argv);

    /*------------------------------------------------------------*/
    Fstk = sf_input("in");
    Fvel = sf_input("velocity");
    Fgam = sf_input("vpvs");
    Fdip = sf_input("dip");
    Fang = sf_output("out");

    if (SF_FLOAT != sf_gettype(Fstk)) sf_error("Need float input");

    az=sf_iaxa(Fstk,1); nz=sf_n(az);
    av=sf_iaxa(Fstk,2); nv=sf_n(av); v0=sf_o(av); dv=sf_d(av);
    ax=sf_iaxa(Fstk,3); nx=sf_n(ax);

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
    dzdx= sf_floatalloc(nz);
    gam = sf_floatalloc(nz);

    sft = fint1_init(ext,nv,0);

    for (ix = 0; ix < nx; ix++) {
	sf_floatread(vel   ,nz   ,Fvel);
	sf_floatread(gam   ,nz   ,Fgam);

	sf_floatread(stk[0],nz*nv,Fstk);
	sf_floatread(dzdx  ,nz   ,Fdip);
	
	/*------------------------------------------------------------*/
	for (iz = 0; iz < nz; iz++) {
	    for (iv = 0; iv < nv; iv++) {
		tmp[iv] = stk[iv][iz];
	    }
	    fint1_set(sft,tmp);

	    v =  vel[iz];
	    g =  gam[iz];
	    d = dzdx[iz];

	    v = v * sqrtf(1+d*d); /* dip correction */
	    e = (1-g)*(1-g) / (4*g);

	    for (ia=0; ia < na; ia++) {
		a = a0+ia*da;          /* angle */
		c = cosf(a/180*SF_PI); /* cos */

		n = v / sqrtf(g) / sqrtf(c*c + e);

		f = (n - v0) / dv;
		fint = f;

		if (fint >= 0 && fint < nv) {
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
