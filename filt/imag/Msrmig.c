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
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */  
    int   nr;             /* number of reference velocities */
    float dt;             /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    bool cw;              /* converted waves flag */

    axa amz,amx,amy;
    axa     alx,aly;
    axa aw,ae;
    axa aj;

    sf_file Fs_s,Fs_r;/*  slowness file S (nlx,nly,nz) */
    sf_file Fw_s,Fw_r;/* wavefield file W ( nx, ny,nw) */
    sf_file Fi;       /*     image file R ( nx, ny,nz) */

    fslice wfl_s,wfl_r,imag;
    fslice slo_s,slo_r;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* converted waves flag */
    if (NULL != sf_getstring("sls")) {
	cw=true;
    } else {
	cw=false;
    }

    if (!sf_getbool("verb",&verb)) verb =  true; /* verbosity flag */
    if (!sf_getfloat("eps",&eps ))  eps =  0.01; /* stability parameter */
    if (!sf_getint(   "nr",&nr  ))   nr =     1; /* maximum number of refs */
    if (!sf_getfloat( "dt",&dt  ))   dt = 0.004; /* time error */
    if (!sf_getint(  "pmx",&pmx ))  pmx =     0; /* padding on x */
    if (!sf_getint(  "pmy",&pmy ))  pmy =     0; /* padding on y */
    if (!sf_getint(  "tmx",&tmx ))  tmx =     0; /* taper on x   */
    if (!sf_getint(  "tmy",&tmy ))  tmy =     0; /* taper on y   */

    /*------------------------------------------------------------*/
    /* SLOWNESS */
    ;      Fs_s = sf_input("slo");
    if(cw) Fs_r = sf_input("sls");
    iaxa(Fs_s,&alx,1); alx.l="lx";
    iaxa(Fs_s,&aly,2); aly.l="ly";
    iaxa(Fs_s,&amz,3); amz.l="mz";
    /* test here if slo and sls have similar sizes */

    ;      slo_s = fslice_init(alx.n*aly.n, amz.n, sizeof(float));
    if(cw) slo_r = fslice_init(alx.n*aly.n, amz.n, sizeof(float));
    ;      fslice_load(Fs_s,slo_s,SF_FLOAT);
    if(cw) fslice_load(Fs_r,slo_r,SF_FLOAT);

    /*------------------------------------------------------------*/    
    /* WAVEFIELD/IMAGE */

    Fw_s = sf_input ( "in");
    Fw_r = sf_input ("rwf");
    Fi   = sf_output("out"); sf_settype(Fi,SF_FLOAT);

    if (SF_COMPLEX != sf_gettype(Fw_s)) sf_error("Need complex   source data");
    if (SF_COMPLEX != sf_gettype(Fw_r)) sf_error("Need complex receiver data");

    aj.n=1; aj.o=0; aj.d=1; aj.l=" ";
    
    iaxa(Fw_s,&amx,1); amx.l="mx"; oaxa(Fi,&amx,1);
    iaxa(Fw_s,&amy,2); amy.l="my"; oaxa(Fi,&amy,2);
    iaxa(Fw_s,&aw,3);  aw.l="w";   oaxa(Fi,&amz,3);
    iaxa(Fw_s,&ae,4);  ae.l="e";   oaxa(Fi,&aj, 4); /* experiments */
    ;                              oaxa(Fi,&aj, 5);

    /* slice management (temp files) */
    wfl_s = fslice_init( amx.n * amy.n, aw.n*ae.n,sizeof(float complex));
    wfl_r = fslice_init( amx.n * amy.n, aw.n*ae.n,sizeof(float complex));
    imag  = fslice_init( amx.n * amy.n, amz.n,    sizeof(float));

    fslice_load(Fw_s,wfl_s,SF_COMPLEX);
    fslice_load(Fw_r,wfl_r,SF_COMPLEX);
    /*------------------------------------------------------------*/
    /* MIGRATION */
    srmig_init (verb,eps,dt,
		ae,aw,amx,amy,amz,alx,aly,
		tmx,tmy,pmx,pmy);

    if(cw) { 
	srmig_cw_init (dt,nr,slo_s,slo_r);
	srmig_cw      (wfl_s,wfl_r,imag);
	srmig_cw_close();
    } else { 
	srmig_pw_init (dt,nr,slo_s);
	srmig_pw      (wfl_s,wfl_r,imag);
	srmig_pw_close();
    }

    srmig_close();

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
