/* 3-D S/R modeling with extended split-step */
/*
  Copyright (C) 2007 Colorado School of Mines
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "srmod3.h"
#include "taper3.h"
#include "ssr3.h"
#include "slow3.h"

#include "weutil.h"
/*^*/

int main (int argc, char *argv[])
{
    bool verb;            /* verbosity */
    bool twoway;          /* two-way traveltime */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    bool cw;              /* converted waves flag */

    sf_axis amz,amx,amy;
    sf_axis alx,aly;
    sf_axis aw,ae;

    int n,nz,nw;

    /* I/O files */
    sf_file Fs_s=NULL,Fs_r=NULL;  /*  slowness file S      (nlx,nly,nz) */
    sf_file Fw_s=NULL,Fw_r=NULL;  /* wavefield file D or U ( nx, ny,nw) */
    sf_file Fr=NULL;              /* reflectivity */

    /* I/O slices */
    sf_fslice wfl_s=NULL,wfl_r=NULL;
    sf_fslice slo_s=NULL,slo_r=NULL;
    sf_fslice refl=NULL;

    int ompchunk=1;
    int ompnth=1;
#ifdef _OPENMP
    int ompath=1; 
#endif

    cub3d cub; /* wavefield hypercube */
    tap3d tap; /* tapering */
    ssr3d ssr; /* SSR operator */
    slo3d s_s; /* slowness */ 
    slo3d s_r; /* slowness */

    ssroperator3d weop;

    float dsmax;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    /* OpenMP data chunk size */
#ifdef _OPENMP
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;
    /* OpenMP available threads */
    
#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    /* converted waves flag */
    if (NULL != sf_getstring("sls")) {
	cw=true;
    } else {
	cw=false;
    }

    if (!sf_getbool("verb",  &verb))    verb =  true; /* verbosity flag */
    if (!sf_getfloat("eps",  &eps ))     eps =  0.01; /* stability parameter */
    if (!sf_getint(  "nrmax",&nrmax))  nrmax =     1; /* maximum number of refs */
    if (!sf_getfloat("dtmax",&dtmax))  dtmax = 0.004; /* time error */
    if (!sf_getint(  "pmx",  &pmx ))     pmx =     0; /* padding on x */
    if (!sf_getint(  "pmy",  &pmy ))     pmy =     0; /* padding on y */
    if (!sf_getint(  "tmx",  &tmx ))     tmx =     0; /* taper on x   */
    if (!sf_getint(  "tmy",  &tmy ))     tmy =     0; /* taper on y   */
    if (!sf_getbool("twoway",&twoway)) twoway= false; /* two-way traveltime */

    /*------------------------------------------------------------*/
    /* SLOWNESS */
    ;      Fs_s = sf_input("slo");
    if(cw) Fs_r = sf_input("sls");
    else   Fs_r = sf_input("slo");

    alx = sf_iaxa(Fs_s,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs_s,2); sf_setlabel(aly,"ly");
    amz = sf_iaxa(Fs_s,3); sf_setlabel(amz,"z" );
    /* test here if slo and sls have similar sizes */

    n  = sf_n(alx)*sf_n(aly);
    nz = sf_n(amz);

    slo_s = sf_fslice_init(n, nz, sizeof(float));
    slo_r = sf_fslice_init(n, nz, sizeof(float));
    sf_fslice_load(Fs_s,slo_s,SF_FLOAT);
    sf_fslice_load(Fs_r,slo_r,SF_FLOAT);

    /*------------------------------------------------------------*/
    /* WAVEFIELD/IMAGE */

    Fw_s = sf_input ( "in");
    Fw_r = sf_output("out"); sf_settype(Fw_r,SF_COMPLEX);
    Fr   = sf_input ("ref");

    if (SF_COMPLEX !=sf_gettype(Fw_s)) sf_error("Need complex source wavefield");
    
    amx = sf_iaxa(Fw_s,1); sf_setlabel(amx,"x"); sf_oaxa(Fw_r,amx,1);
    amy = sf_iaxa(Fw_s,2); sf_setlabel(amy,"y"); sf_oaxa(Fw_r,amy,2);
    aw  = sf_iaxa(Fw_s,3); sf_setlabel(aw, "w"); sf_oaxa(Fw_r,aw,3);
    ae  = sf_iaxa(Fw_s,4); sf_setlabel(ae, "e"); /* experiments */

    n  = sf_n(amx)*sf_n(amy);
    nw = sf_n(aw)*sf_n(ae);

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    wfl_s = sf_fslice_init(n,nw,sizeof(sf_complex));
    wfl_r = sf_fslice_init(n,nw,sizeof(sf_complex));
    refl  = sf_fslice_init(n,nz,sizeof(float));

    sf_fslice_load(Fw_s,wfl_s,SF_COMPLEX);
    sf_fslice_load(Fr,  refl, SF_FLOAT);

    /*------------------------------------------------------------*/
    /* wavefield hypercube */
    cub = srmod3_cube(verb,
		      amx,amy,amz,
		      alx,aly,
		      aw,
		      ae,
		      eps,
		      ompnth,
		      ompchunk);
    
    dsmax = dtmax/cub->amz.d;

    /*------------------------------------------------------------*/
    /* init structures */
    tap = taper_init(cub->amx.n,
		     cub->amy.n,
		     1,
		     SF_MIN(tmx,cub->amx.n-1), /* tmx */
		     SF_MIN(tmy,cub->amy.n-1), /* tmy */
		     0,                        /* tmz */
		     true,true,false);

    ssr = ssr3_init(cub,pmx,pmy,tmx,tmy,dsmax);
    
    s_s = slow3_init(cub,slo_s,nrmax,dsmax,twoway);
    s_r = slow3_init(cub,slo_r,nrmax,dsmax,twoway);

    /*------------------------------------------------------------*/
    /* MODELING */
    weop = srmod3_init(cub);
    
    srmod3(weop,  /* shot-record migration operator */
	   cub,   /* wavefield hypercube dimensions */
	   ssr,   /* SSR operator */
	   tap,   /* tapering operator */
	   s_s,   /* source slowness */
	   s_r,   /* receiver slowness */
	   wfl_s, /* source wavefield */
	   wfl_r, /* receiver wavefield */
	   refl   /* reflectivity */
	);
	   
    srmod3_close(weop);

    /*------------------------------------------------------------*/
    /* close structures   */
    slow3_close(s_s);
    slow3_close(s_r);
    ssr3_close(ssr);
    taper2d_close(tap);

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    if(verb) sf_warning("dump data");
    sf_fslice_dump(Fw_r,wfl_r,SF_COMPLEX);

    /*------------------------------------------------------------*/
    sf_fslice_close(slo_s);
    sf_fslice_close(slo_r);
    sf_fslice_close(wfl_s);
    sf_fslice_close(wfl_r);
    sf_fslice_close(refl);
    
    /*------------------------------------------------------------*/


    exit (0);
}
