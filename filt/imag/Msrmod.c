/* 3-D S/R modeling with extended split-step. */
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
#include "srmod.h"

int main (int argc, char *argv[])
{
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */  
    int   nr;             /* number of reference velocities */
    float dt;             /* time error */
    int   px,py;          /* padding in the k domain */
    int   tx,ty;          /* boundary taper size */

    axa az,ax,ay,aw,alx,aly,ae;

    sf_file Fs;     /*            slowness file S (nlx,nly,nz) */
    sf_file Fd;     /* downgoing wavefield file D ( nx, ny,nw) */
    sf_file Fu;     /*   upgoing wavefield file U ( nx, ny,nw) */
    sf_file Fr;     /* reflectivity */

    fslice slow,dwfl,uwfl,wfld,refl;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (!sf_getbool("verb",&verb)) verb =  true; /* verbosity flag */
    if (!sf_getfloat("eps",&eps ))  eps =  0.01; /* stability parameter */
    if (!sf_getint(   "nr",&nr  ))   nr =     1; /* maximum number of refs */
    if (!sf_getfloat( "dt",&dt  ))   dt = 0.004; /* time error */
    if (!sf_getint(   "px",&px  ))   px =     0; /* padding on x */
    if (!sf_getint(   "py",&py  ))   py =     0; /* padding on y */
    if (!sf_getint(   "tx",&tx  ))   tx =     0; /* taper on x   */
    if (!sf_getint(   "ty",&ty  ))   ty =     0; /* taper on y   */
    
    /* slowness parameters */
    Fs = sf_input ("slo");
    iaxa(Fs,&alx,1); alx.l="lx";
    iaxa(Fs,&aly,2); aly.l="ly";
    iaxa(Fs,&az ,3);  az.l= "z";

    Fd = sf_input ( "in");
    Fu = sf_output("out"); sf_settype(Fu,SF_COMPLEX);
    Fr = sf_input ("ref");

    if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex source wavefield");
    
    iaxa(Fd,&ax,1); ax.l="x"; oaxa(Fu,&ax,1);
    iaxa(Fd,&ay,2); ay.l="y"; oaxa(Fu,&ay,2);
    iaxa(Fd,&aw,3); aw.l="w"; oaxa(Fu,&aw,3);

    /* slice management (temp files) */
    slow = fslice_init(alx.n*aly.n, az.n,sizeof(float));
    dwfl = fslice_init( ax.n* ay.n, aw.n,sizeof(float complex));
    uwfl = fslice_init( ax.n* ay.n, aw.n,sizeof(float complex));
    wfld = fslice_init( ax.n* ay.n, az.n,sizeof(float complex)); /* temp */
    refl = fslice_init( ax.n* ay.n, az.n,sizeof(float));
    fslice_load(Fs,slow,SF_FLOAT);
    fslice_load(Fd,dwfl,SF_COMPLEX);
    fslice_load(Fr,refl,SF_FLOAT);

    srmod_init (verb,eps,dt,
		ae,
		az,aw,
		ax,ay,
		alx,aly,
		tx,ty,
		px,py,
		nr,slow);
    srmod_aloc();
    srmod(dwfl,uwfl,refl,wfld);
    srmod_free();
    srmod_close();
    
    /* slice management (temp files) */
    fslice_dump(Fu,uwfl,SF_COMPLEX);
    fslice_close(slow);
    fslice_close(dwfl);
    fslice_close(uwfl);
    fslice_close(wfld);
    fslice_close(refl);
    
    exit (0);
}
