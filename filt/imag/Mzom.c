/* 3-D zero-offset modeling/migration with extended split-step. */
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
#include "zom.h"

int main (int argc, char *argv[])
{
    char *mode;           /* mode of operation */
    bool inv;             /* forward or adjoint */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant   */  

    int   nr;             /* number of reference velocities */
    float dt;             /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */


    axa az,amx,amy,aw,alx,aly;

    sf_file Fs;    /*  slowness file S(nlx,nly,    nz   ) */
    sf_file Fi;    /*     image file R(nmx,nmy,nhx,nz   ) */
    sf_file Fd,Fe; /*      data file D(nmx,nmy,nhx,   nw) */
    sf_file Fw;    /* wavefield file W(nmx,nmy,nhx,nz,nw) */

    slice slow;
    slice imag;
    slice data,wfld;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* default mode is migration/modeling */
    if (NULL == (mode = sf_getstring("mode"))) mode = "m";

    if (!sf_getbool( "inv",&inv ))  inv = false; /* y=modeling; n=migration */
    if (!sf_getbool("verb",&verb)) verb = false; /* verbosity flag */
    if (!sf_getfloat("eps",&eps ))  eps =  0.01; /* stability parameter */
    if (!sf_getint(   "nr",&nr  ))   nr =     1; /* maximum number of references */
    if (!sf_getfloat( "dt",&dt  ))   dt = 0.004; /* time error */
    if (!sf_getint(  "pmx",&pmx ))  pmx =     0; /* padding on i-line wavenumber */
    if (!sf_getint(  "pmy",&pmy ))  pmy =     0; /* padding on x-line wavenumber */

    if (!sf_getint(  "tmx",&tmx ))  tmx =     0; /* taper size */
    if (!sf_getint(  "tmy",&tmy ))  tmy =     0; /* taper size */

    /* slowness parameters */
    Fs = sf_input ("slowness");
    iaxa(Fs,&alx,1);
    iaxa(Fs,&aly,2);
    iaxa(Fs,&az ,3);
    slow = slice_init(Fs,alx.n,aly.n,az.n);
    
    switch(mode[0]) {
	case 'w': /* save wavefield */
	    Fd = sf_input ( "in");
	    Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
	    if (SF_COMPLEX != sf_gettype(Fd)) sf_error("Need complex data");
 
	    iaxa(Fd,&amx,1); oaxa(Fw,&amx,1);
	    iaxa(Fd,&amy,2); oaxa(Fw,&amy,2);
	    ;                oaxa(Fw,&az ,3);
	    iaxa(Fd,&aw ,3); oaxa(Fw,&aw ,4);

	    data = slice_init(Fd,amx.n,amy.n,     aw.n);
	    wfld = slice_init(Fw,amx.n*amy.n,az.n,aw.n);
	    break;
	case 'd':
	    if (inv) { /*   upward continuation */
		Fe = sf_input ( "in");
		Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
		if (SF_COMPLEX != sf_gettype(Fe)) sf_error("Need complex data");		

		iaxa(Fe,&amx,1); oaxa(Fd,&amx,1);
		iaxa(Fe,&amy,2); oaxa(Fd,&amy,2);
		iaxa(Fe,&aw ,3); oaxa(Fd,&aw ,3);

	    } else {   /* downward continuation */
		Fd = sf_input ( "in");
		Fe = sf_output("out"); sf_settype(Fe,SF_COMPLEX);
		if (SF_COMPLEX != sf_gettype(Fd)) sf_error("Need complex data");
		
		iaxa(Fd,&amx,1); oaxa(Fe,&amx,1);
		iaxa(Fd,&amy,2); oaxa(Fe,&amy,2);
		iaxa(Fd,&aw ,3); oaxa(Fe,&aw ,3);
	    }
	    data = slice_init(Fd,amx.n,amy.n,aw.n);
	    wfld = slice_init(Fe,amx.n,amy.n,aw.n);
	    break;
	case 'm':
	default:
	    if (inv) { /* modeling */
		Fi = sf_input ( "in");
		Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
		if (SF_FLOAT != sf_gettype(Fi)) sf_error("Need float image");
		
		if (!sf_getint  ("nw",&aw.n)) sf_error ("Need nw=");
		if (!sf_getfloat("dw",&aw.d)) sf_error ("Need dw=");
		if (!sf_getfloat("w0",&aw.o)) aw.o=0.;
		
		iaxa(Fi,&amx,1); oaxa(Fd,&amx,1);
		iaxa(Fi,&amy,2); oaxa(Fd,&amy,2);
		iaxa(Fi,&az ,3); oaxa(Fd,&aw ,3);
		
	    } else { /* migration */
		Fd = sf_input ( "in");
		Fi = sf_output("out"); sf_settype(Fi,SF_FLOAT);
		if (SF_COMPLEX != sf_gettype(Fd)) sf_error("Need complex data");
		
		iaxa(Fd,&amx,1); oaxa(Fi,&amx,1);
		iaxa(Fd,&amy,2); oaxa(Fi,&amy,2);
		iaxa(Fd,&aw ,3); oaxa(Fi,&az ,3);
	    }
	    data = slice_init(Fd,amx.n,amy.n,aw.n);
	    imag = slice_init(Fi,amx.n,amy.n,az.n);
	    break;
    }
    
    zom_init (verb,eps,dt,
	      az,aw,
	      amx,amy,
	      alx,aly,
	      tmx,tmy,
	      pmx,pmy,
	      nr,slow);

    switch(mode[0]) {
	case 'w':
	    cawfl(    data,wfld);
	    break;
	case 'd':
	    cadtm(inv,data,wfld);
	    break;
	case 'm':
	default:
	    zomig_init();
	    zomig(inv,data,imag);
	    zomig_close();
	    break;
    }

    zom_close();
    
    exit (0);
}
