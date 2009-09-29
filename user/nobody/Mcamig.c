/* 3-D common-azimuth modeling/migration with extended split-step */
/*
  Copyright (C) 2006 Colorado School of Mines
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
#include "camig.h"

int main (int argc, char *argv[])
{
    char *mode;           /* mode of operation */
    bool inv;             /* forward or adjoint */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy,phx;    /* padding in the k domain */
    int   tmx,tmy,thx;    /* boundary taper size */

    sf_axis az,amx,amy,aw,alx,aly,ahx,ae;
    int nz, n, nw;
    float dw, w0;

    sf_file Fs;           /*  slowness file S(nlx,nly,    nz   ) */
    sf_file Fi=NULL;      /*     image file R(nmx,nmy,nhx,nz   ) */
    sf_file Fd,Fw=NULL;   /*      data file D(nmx,nmy,nhx,   nw) */

    fslice slow,imag=NULL,data,wfld=NULL;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* default mode is migration/modeling */
    if (NULL == (mode = sf_getstring("mode"))) mode = "m";

    if (!sf_getbool( "inv",  &inv  ))   inv = false; /* y=modeling; n=migration */
    if (!sf_getbool( "verb", &verb ))  verb = false; /* verbosity flag */
    if (!sf_getfloat("eps",  &eps  ))   eps =  0.01; /* stability parameter */
    if (!sf_getint(  "nrmax",&nrmax)) nrmax =     1; /* maximum number of refs */
    if (!sf_getfloat("dtmax",&dtmax)) dtmax = 0.004; /* time error */
    if (!sf_getint(  "pmx",  &pmx  ))   pmx =     0; /* padding mx*/
    if (!sf_getint(  "pmy",  &pmy  ))   pmy =     0; /* padding my*/
    if (!sf_getint(  "phx",  &phx  ))   phx =     0; /* padding hx*/

    if (!sf_getint(  "tmx",  &tmx  ))   tmx =     0; /* taper mx */
    if (!sf_getint(  "tmy",  &tmy  ))   tmy =     0; /* taper my */
    if (!sf_getint(  "thx",  &thx  ))   thx =     0; /* taper hx */

    /* slowness parameters */
    Fs = sf_input ("slo");
    alx = sf_iaxa(Fs,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs,2); sf_setlabel(aly,"ly");
    az  = sf_iaxa(Fs,3); nz  = sf_n(az);  sf_setlabel(az,"z");
    slow = fslice_init(sf_n(alx)*sf_n(aly),nz,sizeof(float));
    fslice_load(Fs,slow,SF_FLOAT);
    
      switch(mode[0]) {
	case 'w': /* save wavefield */
	    Fd = sf_input ( "in");
	    Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
	    if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
 
	    amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fw,amx,1);
	    amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fw,amy,2);
	    ahx = sf_iaxa(Fd,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fw,ahx,3); 
	    ;                                           sf_oaxa(Fw,az ,4);
	    aw  = sf_iaxa(Fd,4); sf_setlabel(aw ,"w" ); sf_oaxa(Fw,aw ,5);
	    ae = sf_maxa(1,0,1);

	    n = sf_n(amx)*sf_n(amy)*sf_n(ahx);
	    nw = sf_n(aw);

	    data = fslice_init(n,nw,sizeof(sf_complex));
	    wfld = fslice_init(n, nz*nw,sizeof(sf_complex));

	    fslice_load(Fd,data,SF_COMPLEX);

	    break;
	case 'd':
	    if (inv) { /*   upward continuation */
		Fw = sf_input ( "in");
		Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
		if (SF_COMPLEX !=sf_gettype(Fw)) sf_error("Need complex data");

		amx = sf_iaxa(Fw,1); sf_setlabel(amx,"mx"); sf_oaxa(Fd,amx,1);
		amy = sf_iaxa(Fw,2); sf_setlabel(amy,"my"); sf_oaxa(Fd,amy,2);
		ahx = sf_iaxa(Fw,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fd,ahx,3);
		aw  = sf_iaxa(Fw,4); sf_setlabel(aw , "w"); sf_oaxa(Fd,aw ,4);
		ae  = sf_iaxa(Fw,5); sf_setlabel(ae,  "e"); sf_oaxa(Fd,ae ,5);

		n = sf_n(amx)*sf_n(amy)*sf_n(ahx);
		nw = sf_n(aw)*sf_n(ae);

		data = fslice_init(n,nw,sizeof(sf_complex));
		wfld = fslice_init(n,nw,sizeof(sf_complex));

		fslice_load(Fw,wfld,SF_COMPLEX);
	    } else {   /* downward continuation */
		Fd = sf_input ( "in");
		Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
		if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
		
		amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fw,amx,1);
		amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fw,amy,2);
		ahx = sf_iaxa(Fd,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fw,ahx,3);
		aw  = sf_iaxa(Fd,4); sf_setlabel(aw , "w"); sf_oaxa(Fw,aw ,4);
		ae  = sf_iaxa(Fd,5); sf_setlabel(ae , "e"); sf_oaxa(Fw,ae ,5);

		n = sf_n(amx)*sf_n(amy)*sf_n(ahx);
		nw = sf_n(aw)*sf_n(ae);

		data = fslice_init(n,nw,sizeof(sf_complex));
		wfld = fslice_init(n,nw,sizeof(sf_complex));

		fslice_load(Fd,data,SF_COMPLEX);
	    }
	    break;
	case 'm':
	default:
	    if (inv) { /* modeling */
		Fi = sf_input ( "in");
		Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
		if (SF_FLOAT !=sf_gettype(Fi)) sf_error("Need float image");
		
		if (!sf_getint  ("nw",&nw)) sf_error ("Need nw=");
		if (!sf_getfloat("dw",&dw)) sf_error ("Need dw=");
		if (!sf_getfloat("ow",&w0)) w0=0.;
		aw = sf_maxa(nw,w0,dw); 
		sf_setlabel(aw, "w");
		sf_setunit(aw, "1/s");
		
		amx = sf_iaxa(Fi,1); sf_setlabel(amx,"mx"); sf_oaxa(Fd,amx,1);
		amy = sf_iaxa(Fi,2); sf_setlabel(amy,"my"); sf_oaxa(Fd,amy,2);
		ahx = sf_iaxa(Fi,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fd,ahx,3);
		az  = sf_iaxa(Fi,4); sf_setlabel(az , "z"); sf_oaxa(Fd,aw ,4);
		ae = sf_maxa(1,0,1);
		
		n = sf_n(amx)*sf_n(amy)*sf_n(ahx);

		data = fslice_init(n,nw,sizeof(sf_complex));
		imag = fslice_init(n,nz,sizeof(float));

		fslice_load(Fi,imag,SF_FLOAT);
	    } else { /* migration */
		Fd = sf_input ( "in");
		Fi = sf_output("out"); sf_settype(Fi,SF_FLOAT);
		if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
		
		amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fi,amx,1);
		amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fi,amy,2);
		ahx = sf_iaxa(Fd,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fi,ahx,3);
		aw  = sf_iaxa(Fd,4); sf_setlabel(aw , "w"); sf_oaxa(Fi,az ,4);
		ae = sf_maxa(1,0,1);

		n = sf_n(amx)*sf_n(amy)*sf_n(ahx);
		nw = sf_n(aw);

		data = fslice_init(n,nw,sizeof(sf_complex));
		imag = fslice_init(n,nz,sizeof(float));
	    
		fslice_load(Fd,data,SF_COMPLEX);
	    } 
	    break;
    }

    /*------------------------------------------------------------*/
    
    camig_init(verb,eps,dtmax,
	       az,aw,ae,
	       amx,amy,ahx,
	       alx,aly,
	       tmx,tmy,thx,
	       pmx,pmy,phx,
	       nrmax,slow);
    
    switch(mode[0]) {
	case 'w':
	    cawfl(    data,wfld);
	    break;
	case 'd':
	    cadtm(inv,data,wfld);
	    break;
	case 'm':
	default:
	    camig_aloc();
	    camig(inv,data,imag);
	    camig_free();
	    break;
    }

    camig_close();

    /*------------------------------------------------------------*/

    switch(mode[0]) {
	case 'w':
	    fslice_dump(Fw,wfld,SF_COMPLEX);
	    fslice_close(data);
	    fslice_close(wfld);
	    break;
	case 'd':
	    if(inv) fslice_dump(Fd,data,SF_COMPLEX);
	    else    fslice_dump(Fw,wfld,SF_COMPLEX);
	    fslice_close(data);
	    fslice_close(wfld);
	    break;
	case 'm':
	    if(inv) fslice_dump(Fd,data,SF_COMPLEX);
	    else    fslice_dump(Fi,imag,SF_FLOAT);
	    fslice_close(data);
	    fslice_close(imag);
	default:
	    break;
    }
    fslice_close(slow);

    exit (0);
}
