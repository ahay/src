/* 3-D S/R migration with extended split-step. */
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
#include "srmig.h"

int main (int argc, char *argv[])
{
    bool inv;             /* forward or adjoint */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */  
    int   nr;             /* number of reference velocities */
    float dt;             /* time error */
    int   px,py;          /* padding in the k domain */
    int   tx,ty;          /* boundary taper size */

    axa az,ax,ay,aw,alx,aly;
    axa ae;
    axa aj;

    sf_file Fs;     /*           slowness file S (nlx,nly,nz) */
    sf_file Fi;     /*              image file R ( nx, ny,nz) */
    sf_file Fus;    /*   source wavefield file Us( nx, ny,nw) */
    sf_file Fur;    /* receiver wavefield file Us( nx, ny,nw) */

    slice slow;
    slice imag;
    slice sdat,rdat;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (!sf_getbool( "inv",&inv ))  inv = false; /* y=modeling; n=migration */
    if (!sf_getbool("verb",&verb)) verb =  true; /* verbosity flag */
    if (!sf_getfloat("eps",&eps ))  eps =  0.01; /* stability parameter */
    if (!sf_getint(   "nr",&nr  ))   nr =     1; /* maximum number of references */
    if (!sf_getfloat( "dt",&dt  ))   dt = 0.004; /* time error */
    if (!sf_getint(   "px",&px  ))   px =     0; /* padding on i-line wavenumber */
    if (!sf_getint(   "py",&py  ))   py =     0; /* padding on x-line wavenumber */
    if (!sf_getint(   "tx",&tx  ))   tx =     0; /* taper size */
    if (!sf_getint(   "ty",&ty  ))   ty =     0; /* taper size */

    /* slowness parameters */
    Fs = sf_input ("slo");
    iaxa(Fs,&alx,1); alx.l="lx";
    iaxa(Fs,&aly,2); aly.l="ly";
    iaxa(Fs,&az ,3);  az.l= "z";
    slow = slice_init(Fs,alx.n,aly.n,az.n);
    
    Fus = sf_input ( "in");
    Fur = sf_input ("rwf");
    Fi  = sf_output("out"); sf_settype(Fi,SF_FLOAT);

    if (SF_COMPLEX != sf_gettype(Fus)) sf_error("Need complex   source data");
    if (SF_COMPLEX != sf_gettype(Fur)) sf_error("Need complex receiver data");

    aj.n=1; aj.o=0; aj.d=1; aj.l=" ";
    
    iaxa(Fus,&ax,1); ax.l="x"; oaxa(Fi,&ax,1);
    iaxa(Fus,&ay,2); ay.l="y"; oaxa(Fi,&ay,2);
    iaxa(Fus,&aw,3); aw.l="w"; oaxa(Fi,&az,3);
    iaxa(Fus,&ae,4); ae.l="e"; oaxa(Fi,&aj,4); /* no of experiments */
    ;                          oaxa(Fi,&aj,5);

    sdat = slice_init(Fus,ax.n,ay.n,aw.n);
    rdat = slice_init(Fur,ax.n,ay.n,aw.n);
    imag = slice_init( Fi,ax.n,ay.n,az.n);

    srmig_init (verb,eps,dt,
		ae,
		az,aw,
		ax,ay,
		alx,aly,
		tx,ty,
		px,py,
		nr,slow);

    srmig_aloc();

    srmig(inv,sdat,rdat,imag);

    srmig_free();

    srmig_close();
	
    exit (0);
}
