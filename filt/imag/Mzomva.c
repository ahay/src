/* 3-D zero-offset MVA */
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
#include "zomva.h"

int main (int argc, char *argv[])
{
    
    bool inv;             /* forward or adjoint */
    bool twoway;          /* two-way traveltime */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    
    axa az,amx,amy,aw;
    axa    alx,aly;

    sf_file Bs;  /*  slowness  file S (nlx,nly,nz)    */
    sf_file Bw;  /*  wavefield file W (nmx,nmy,nz,nw) */
    sf_file Ps;  /* slowness perturbation dS (nmx,nmy,nz) complex */
    sf_file Pi;  /*    image perturbation dI (nmx,nmy,nz) complex */

    fslice Bslow,Bwfld;
    fslice Pslow,Pimag;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (!sf_getbool(  "verb",&verb ))   verb = false; /* verbosity flag */
    if (!sf_getfloat(  "eps",&eps  ))    eps =  0.01; /* stability parameter */
    if (!sf_getbool(   "inv",&inv  ))    inv = false; /* y=modeling; n=migration */
    if (!sf_getbool("twoway",&twoway)) twoway=  true; /* two-way traveltime */
    if (!sf_getint(  "nrmax",&nrmax))  nrmax =     1; /* maximum number of references */
    if (!sf_getfloat("dtmax",&dtmax))  dtmax = 0.004; /* time error */
    if (!sf_getint(    "pmx",&pmx  ))    pmx =     0; /* padding on x */
    if (!sf_getint(    "pmy",&pmy  ))    pmy =     0; /* padding on y*/
    if (!sf_getint(    "tmx",&tmx  ))    tmx =     0; /* taper on x*/
    if (!sf_getint(    "tmy",&tmy  ))    tmy =     0; /* taper on y */

    /* slowness parameters */
    Bs = sf_input ("slo");
    iaxa(Bs,&alx,1); alx.l="lx";
    iaxa(Bs,&aly,2); aly.l="ly";
    iaxa(Bs,&az ,3);  az.l= "z";
    Bslow = fslice_init(alx.n*aly.n,az.n,sizeof(float));
    fslice_load(Bs,Bslow,SF_FLOAT);

    /* wavefield parameters */
    Bw = sf_input("wfl");
    if (SF_COMPLEX !=sf_gettype(Bw)) sf_error("Need complex wavefield");
    
    iaxa(Bw,&amx,1); amx.l="mx";
    iaxa(Bw,&amy,2); amy.l="my";
    iaxa(Bw,&az, 3);  az.l= "z";
    iaxa(Bw,&aw ,4);  aw.l= "w";

    Bwfld = fslice_init(amx.n*amy.n,az.n*aw.n,sizeof(float complex));
    fslice_load(Bw,Bwfld,SF_COMPLEX);

    if (inv) { /* adjoint: image -> slowness */
	
	Pi = sf_input ( "in");
	if (SF_COMPLEX !=sf_gettype(Pi)) sf_error("Need complex image perturbation");

	Ps = sf_output("out"); sf_settype(Ps,SF_COMPLEX);
	oaxa(Ps,&amx,1);
	oaxa(Ps,&amy,2);
	oaxa(Ps,&az, 3);

	Pslow = fslice_init(amx.n*amy.n, az.n, sizeof(float complex));

	Pimag = fslice_init(amx.n*amy.n, az.n, sizeof(float complex));
	fslice_load(Pi,Pimag,SF_COMPLEX);

    } else {   /* forward: slowness -> image */

	Ps = sf_input ( "in");
	if (SF_COMPLEX !=sf_gettype(Ps)) sf_error("Need complex slowness perturbation");

	Pi = sf_output("out"); sf_settype(Pi,SF_COMPLEX);
	oaxa(Pi,&amx,1);
	oaxa(Pi,&amy,2);
	oaxa(Pi,&az, 3);
	
	Pslow = fslice_init(amx.n*amy.n, az.n, sizeof(float complex));
	fslice_load(Ps,Pslow,SF_COMPLEX);

	Pimag = fslice_init(amx.n*amy.n, az.n, sizeof(float complex));

    }
    /*------------------------------------------------------------*/

    zomva_init(verb,eps,twoway,dtmax,
	       az,aw,
	       amx,amy,
	       alx,aly,
	       tmx,tmy,
	       pmx,pmy,
	       nrmax,Bslow,Bwfld);

    zomva_aloc();
    zomva(inv,Pslow,Pimag);
    zomva_free();

    zomva_close();


    /*------------------------------------------------------------*/
    if(inv) fslice_dump(Ps,Pslow,SF_COMPLEX);
    else    fslice_dump(Pi,Pimag,SF_COMPLEX);
    fslice_close(Pimag);
    fslice_close(Pslow);

    fslice_close(Bwfld);
    fslice_close(Bslow);

    exit (0);
}
