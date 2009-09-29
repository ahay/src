/* 3-D S/R modeling with extended split-step */
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
#include "srmod.h"

int main (int argc, char *argv[])
{
    bool verb;            /* verbosity */
    bool incore;          /* in core execution */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    bool cw;              /* converted waves flag */

    sf_axis az,ax,ay,aw,alx,aly;

    int n,nz,nw;

    sf_file Fs_s=NULL,Fs_r=NULL;  /*  slowness file S      (nlx,nly,nz) */
    sf_file Fw_s=NULL,Fw_r=NULL;  /* wavefield file D or U ( nx, ny,nw) */
    sf_file Fr=NULL;              /* reflectivity */

    fslice wfl_s=NULL,wfl_r=NULL,refl=NULL;
    fslice slo_s=NULL,slo_r=NULL;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* converted waves flag */
    if (NULL != sf_getstring("sls")) {
	cw=true;
    } else {
	cw=false;
    }

    if (!sf_getbool("verb",  &verb))    verb =  true; /* verbosity flag */
    if (!sf_getbool("incore",&incore))incore = false; /* in core execution */
    if (!sf_getfloat("eps",  &eps ))     eps =  0.01; /* stability parameter */
    if (!sf_getint(  "nrmax",&nrmax))  nrmax =     1; /* maximum number of refs */
    if (!sf_getfloat("dtmax",&dtmax))  dtmax = 0.004; /* time error */
    if (!sf_getint(  "pmx",  &pmx ))     pmx =     0; /* padding on x */
    if (!sf_getint(  "pmy",  &pmy ))     pmy =     0; /* padding on y */
    if (!sf_getint(  "tmx",  &tmx ))     tmx =     0; /* taper on x   */
    if (!sf_getint(  "tmy",  &tmy ))     tmy =     0; /* taper on y   */
    
    /*------------------------------------------------------------*/
    /* SLOWNESS */
    ;      Fs_s = sf_input("slo");
    if(cw) Fs_r = sf_input("sls");
    alx = sf_iaxa(Fs_s,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs_s,2); sf_setlabel(aly,"ly");
    az  = sf_iaxa(Fs_s,3); sf_setlabel(az , "z");
    /* test here if slo and sls have similar sizes */

    n = sf_n(alx)*sf_n(aly);
    nz = sf_n(az);

    ;      slo_s = fslice_init(n,nz,sizeof(float));
    if(cw) slo_r = fslice_init(n,nz,sizeof(float));
    ;      fslice_load(Fs_s,slo_s,SF_FLOAT);
    if(cw) fslice_load(Fs_r,slo_r,SF_FLOAT);
    
    /*------------------------------------------------------------*/
    /* WAVEFIELD/IMAGE */
    Fw_s = sf_input ( "in");
    Fw_r = sf_output("out"); sf_settype(Fw_r,SF_COMPLEX);
    Fr   = sf_input ("ref");

    if (SF_COMPLEX !=sf_gettype(Fw_s)) 
	sf_error("Need complex source wavefield");
    
    ax = sf_iaxa(Fw_s,1); sf_setlabel(ax,"x"); sf_oaxa(Fw_r,ax,1);
    ay = sf_iaxa(Fw_s,2); sf_setlabel(ay,"y"); sf_oaxa(Fw_r,ay,2);
    aw = sf_iaxa(Fw_s,3); sf_setlabel(aw,"w"); sf_oaxa(Fw_r,aw,3);

    n = sf_n(ax)*sf_n(ay);
    nw = sf_n(aw);

    /* slice management (temp files) */
    wfl_s = fslice_init(n,nw,sizeof(sf_complex));
    wfl_r = fslice_init(n,nw,sizeof(sf_complex));
    refl  = fslice_init(n,nz,sizeof(float));

    fslice_load(Fw_s,wfl_s,SF_COMPLEX);
    fslice_load(Fr,  refl, SF_FLOAT);

    /*------------------------------------------------------------*/
    /* MODELING */
    srmod_init (verb,incore,eps,dtmax,
		az,aw,ax,ay,alx,aly,
		tmx,tmy,pmx,pmy);
    
    if(cw) { 
	srmod_cw_init (dtmax,nrmax,slo_s,slo_r);
	srmod_cw      (wfl_s,wfl_r,refl);
	srmod_cw_close();
    } else { 
	srmod_pw_init (dtmax,nrmax,slo_s);
	srmod_pw      (wfl_s,wfl_r,refl);
	srmod_pw_close();
    }

    srmod_close();

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    fslice_dump(Fw_r,wfl_r,SF_COMPLEX);
    ;      fslice_close(slo_s);
    if(cw) fslice_close(slo_r);
    ;      fslice_close(wfl_s);
    ;      fslice_close(wfl_r);
    ;      fslice_close(refl);

    exit (0);
}
