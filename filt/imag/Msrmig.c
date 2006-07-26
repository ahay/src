/* 3-D S/R migration with extended split-step */

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
#include "srmig.h"

#include "img.h"

int main (int argc, char *argv[])
{
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    bool cw;              /* converted waves flag */
    char *itype;          /* imaging type 
			     o = zero offset (default)
			     x = space offset
			     t = time offset
			  */
    bool hsym;

    sf_axis amx,amy,amz;
    sf_axis alx,aly;
    sf_axis aw,ae;
    sf_axis ahx,ahy,ahz,aht;

    int n,nm,nz,nw,nh,nhx,nhy,nhz;
    float oh,dh;

    sf_file Fs_s=NULL,Fs_r=NULL;/*  slowness file S (nlx,nly,nz) */
    sf_file Fw_s=NULL,Fw_r=NULL;/* wavefield file W ( nx, ny,nw) */
    sf_file Fi=NULL;            /*     image file R ( nx, ny,nz) */

    fslice wfl_s=NULL,wfl_r=NULL,imag=NULL;
    fslice slo_s=NULL,slo_r=NULL;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (NULL == (itype = sf_getstring("itype"))) itype = "o";

    /* converted waves flag */
    if (NULL != sf_getstring("sls")) {
	cw=true;
    } else {
	cw=false;
    }

    if (!sf_getbool(  "verb",&verb ))  verb =  true; /* verbosity flag */
    if (!sf_getfloat(  "eps",&eps  ))   eps =  0.01; /* stability parameter */
    if (!sf_getint(  "nrmax",&nrmax)) nrmax =     1; /* max number of refs */
    if (!sf_getfloat("dtmax",&dtmax)) dtmax = 0.004; /* max time error */
    if (!sf_getint(    "pmx",&pmx  ))   pmx =     0; /* padding on x */
    if (!sf_getint(    "pmy",&pmy  ))   pmy =     0; /* padding on y */
    if (!sf_getint(    "tmx",&tmx  ))   tmx =     0; /* taper on x   */
    if (!sf_getint(    "tmy",&tmy  ))   tmy =     0; /* taper on y   */

    /*------------------------------------------------------------*/
    /* SLOWNESS */
    ;      Fs_s = sf_input("slo");
    if(cw) Fs_r = sf_input("sls");
    alx = sf_iaxa(Fs_s,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs_s,2); sf_setlabel(aly,"ly");
    amz = sf_iaxa(Fs_s,3); sf_setlabel(amz,"mz");
    /* test here if slo and sls have similar sizes */

    n = sf_n(alx)*sf_n(aly);
    nz = sf_n(amz);

    ;      slo_s = fslice_init(n, nz, sizeof(float));
    if(cw) slo_r = fslice_init(n, nz, sizeof(float));
    ;      fslice_load(Fs_s,slo_s,SF_FLOAT);
    if(cw) fslice_load(Fs_r,slo_r,SF_FLOAT);

    /*------------------------------------------------------------*/    
    /* WAVEFIELD/IMAGE */

    Fw_s = sf_input ( "in");
    Fw_r = sf_input ("rwf");
    Fi   = sf_output("out"); sf_settype(Fi,SF_FLOAT);

    if (SF_COMPLEX != sf_gettype(Fw_s)) sf_error("Need complex   source data");
    if (SF_COMPLEX != sf_gettype(Fw_r)) sf_error("Need complex receiver data");
    
    amx = sf_iaxa(Fw_s,1); sf_setlabel(amx,"mx"); sf_oaxa(Fi,amx,1);
    amy = sf_iaxa(Fw_s,2); sf_setlabel(amy,"my"); sf_oaxa(Fi,amy,2);
    aw  = sf_iaxa(Fw_s,3); sf_setlabel(aw ,"w" ); sf_oaxa(Fi,amz,3);
    ae  = sf_iaxa(Fw_s,4); sf_setlabel(ae ,"e" );   /* experiments */

    nm = sf_n(amx)*sf_n(amy);
    n = nm*nz;

    switch(itype[0]) {
	case 't': /* time offset imaging condition */
	    if(!sf_getint  ("nht",&nh)) nh=1;
	    if(!sf_getfloat("oht",&oh)) oh=0;
	    if(!sf_getfloat("dht",&dh)) dh=0.1;
	    aht = sf_maxa(nh,oh,dh); sf_setlabel(aht,"ht");
	    
	    sf_oaxa(Fi,aht,4);
	    sf_putint(Fi,"n5",1);

	    imag = fslice_init(n,nh,sizeof(float));
	    imgt_init(amz,amx,amy,aht,aw,imag);

	    break;
	case 'x': /* space offset imaging condition */
	    if(!sf_getint("nhx",&nhx)) nhx=1;
	    if(!sf_getint("nhy",&nhy)) nhy=1;
	    if(!sf_getint("nhz",&nhz)) nhz=1;
	    ahx = sf_maxa(nhx,0.,sf_d(amx)); sf_setlabel(ahx,"hx");
	    ahy = sf_maxa(nhy,0.,sf_d(amy)); sf_setlabel(ahy,"hy");
	    ahz = sf_maxa(nhz,0.,sf_d(amz)); sf_setlabel(ahz,"hz");

	    if(!sf_getbool("hsym",&hsym)) hsym = false;
	    if(hsym) {
		if(nhx>1) {
		    sf_seto(ahx,- nhx*sf_d(ahx)); nhx *= 2; sf_setn(ahx,nhx); }
		if(nhy>1) {
		    sf_seto(ahy,- nhy*sf_d(ahy)); nhx *= 2; sf_setn(ahy,nhy); }
		if(nhz>1) {
		    sf_seto(ahz,- nhz*sf_d(ahz)); nhz *= 2; sf_setn(ahz,nhz); }
	    }

	    sf_oaxa(Fi,ahx,4);
	    sf_oaxa(Fi,ahy,5);
	    sf_oaxa(Fi,ahz,6);

	    imag = fslice_init(n,nhx*nhy*nhz,sizeof(float));
	    imgx_init(amz,amx,amy,ahx,ahy,ahz,imag);

	    break;
	case 'o': /* zero offset imaging condition */
	default:
	    sf_putint(Fi,"n4",1);
	    sf_putint(Fi,"n5",1);

	    imag = fslice_init(n,1,sizeof(float));
	    imgo_init(amz,amx,amy,imag);

	    break;
    }

    /* slice management (temp files) */
    nw = sf_n(aw)*sf_n(ae);

    wfl_s = fslice_init( nm, nw, sizeof(sf_complex));
    wfl_r = fslice_init( nm, nw, sizeof(sf_complex));

    fslice_load(Fw_s,wfl_s,SF_COMPLEX);
    fslice_load(Fw_r,wfl_r,SF_COMPLEX);
    /*------------------------------------------------------------*/
    /* MIGRATION */
    srmig_init (verb,eps,dtmax,
		ae,aw,amx,amy,amz,alx,aly,
		tmx,tmy,pmx,pmy);
    
    if(cw) { 
	srmig_cw_init (dtmax,nrmax,slo_s,slo_r);
	switch(itype[0]) {
	    case 't':          srmig_cw(wfl_s,wfl_r,imag, &imgt); break;
	    case 'x':          srmig_cw(wfl_s,wfl_r,imag, &imgx); break;
	    case 'o': default: srmig_cw(wfl_s,wfl_r,imag, &imgo); break;
	}
	srmig_cw_close();
    } else { 
	srmig_pw_init (dtmax,nrmax,slo_s);
	switch(itype[0]) {
	    case 't':          srmig_pw(wfl_s,wfl_r,imag, &imgt); break;
	    case 'x':          srmig_pw(wfl_s,wfl_r,imag, &imgx); break;
	    case 'o': default: srmig_pw(wfl_s,wfl_r,imag, &imgo); break;
	}
	srmig_pw_close();
    }
    
    srmig_close();

    switch(itype[0]) {
	case 't':          imgt_close(); break;
	case 'x':          imgx_close(); break;
	case 'o': default: imgo_close(); break;
    }

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    fslice_dump(Fi,imag,SF_FLOAT);

    ;      fslice_close(slo_s);
    if(cw) fslice_close(slo_r);
    ;      fslice_close(wfl_s);
    ;      fslice_close(wfl_r);
    ;      fslice_close(imag);

    exit (0);
}
