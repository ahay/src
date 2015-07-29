/* 3-D common-azimuth modeling/migration with extended SSF */

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

#include "camig3.h"
#include "taper3.h"
#include "cam3.h"
#include "slow3.h"

#include "weutil.h"
/*^*/

int main (int argc, char *argv[])
{
    const char *mode;     /* mode of operation */
    bool verb;            /* verbosity */
    bool inv;             /* forward or adjoint */
    bool twoway;          /* two-way traveltime */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy,phx;    /* padding in the k domain */
    int   tmx,tmy,thx;    /* boundary taper size */

    sf_axis amx,amy,amz;
    sf_axis ahx;
    sf_axis alx,aly;
    sf_axis aw;

    int nz, n, nw;
    float dw, w0;

    /* I/O files */
    sf_file Fs=NULL;      /*  slowness file S(nlx,nly,    nz   ) */
    sf_file Fi=NULL;      /*     image file R(nmx,nmy,nhx,nz   ) */
    sf_file Fd=NULL;      /*      data file D(nmx,nmy,nhx,   nw) */
    sf_file Fw=NULL; 

    /* I/O slices */
    sf_fslice slow=NULL;
    sf_fslice imag=NULL;
    sf_fslice data=NULL;
    sf_fslice wfld=NULL;

    int ompchunk=1;
    int ompnth=1;
#ifdef _OPENMP
    int ompath=1; 
#endif
    
    cub3d cub; /* wavefield hypercube */
    tap3d tap; /* tapering */
    cam3d cam; /* CAM operator */
    slo3d slo; /* slowness */

    camoperator3d weop;

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

    /* default mode is migration/modeling */
    if (NULL == (mode = sf_getstring("mode"))) mode = "m";

    if (!sf_getbool( "verb", &verb ))  verb = false; /* verbosity flag */
    if (!sf_getfloat("eps",  &eps  ))   eps =  0.01; /* stability parameter */
    if (!sf_getbool( "inv",  &inv  ))   inv = false; /* y=modeling; n=migration */

    if (!sf_getbool("twoway",&twoway))twoway= false; /* two-way traveltime */
    if (!sf_getint(  "nrmax",&nrmax)) nrmax =     1; /* maximum number of refs */
    if (!sf_getfloat("dtmax",&dtmax)) dtmax = 0.004; /* time error */

    if (!sf_getint(  "pmx",  &pmx  ))   pmx =     0; /* padding mx*/
    if (!sf_getint(  "pmy",  &pmy  ))   pmy =     0; /* padding my*/
    if (!sf_getint(  "phx",  &phx  ))   phx =     0; /* padding hx*/

    if (!sf_getint(  "tmx",  &tmx  ))   tmx =     0; /* taper mx */
    if (!sf_getint(  "tmy",  &tmy  ))   tmy =     0; /* taper my */
    if (!sf_getint(  "thx",  &thx  ))   thx =     0; /* taper hx */

    /* slowness parameters */
    Fs  = sf_input ("slo");
    alx = sf_iaxa(Fs,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs,2); sf_setlabel(aly,"ly");
    amz = sf_iaxa(Fs,3); sf_setlabel(amz,"z" );

    n   = sf_n(alx)*sf_n(aly);
    nz  = sf_n(amz); 

    slow = sf_fslice_init(n,nz,sizeof(float));
    sf_fslice_load(Fs,slow,SF_FLOAT);
    
    switch(mode[0]) {
	case 'w': /* save wavefield */
	    Fd = sf_input ( "in");
	    Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
	    if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
	    
	    amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fw,amx,1);
	    amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fw,amy,2);
	    ahx = sf_iaxa(Fd,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fw,ahx,3); 
	    ;                                           sf_oaxa(Fw,amz,4);
	    aw  = sf_iaxa(Fd,4); sf_setlabel(aw ,"w" ); sf_oaxa(Fw,aw ,5);
	    
	    n  = sf_n(amx)*sf_n(amy)*sf_n(ahx);
	    nw = sf_n(aw);
	    
	    data = sf_fslice_init(n,   nw,sizeof(sf_complex));
	    wfld = sf_fslice_init(n,nz*nw,sizeof(sf_complex));
	    
	    sf_fslice_load(Fd,data,SF_COMPLEX);
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
		
		n  = sf_n(amx)*sf_n(amy)*sf_n(ahx);
		nw = sf_n(aw);
		
		data = sf_fslice_init(n,nw,sizeof(sf_complex));
		wfld = sf_fslice_init(n,nw,sizeof(sf_complex));
		
		sf_fslice_load(Fw,wfld,SF_COMPLEX);
	    } else {   /* downward continuation */
		Fd = sf_input ( "in");
		Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
		if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
		
		amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fw,amx,1);
		amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fw,amy,2);
		ahx = sf_iaxa(Fd,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fw,ahx,3);
		aw  = sf_iaxa(Fd,4); sf_setlabel(aw , "w"); sf_oaxa(Fw,aw ,4);
		
		n  = sf_n(amx)*sf_n(amy)*sf_n(ahx);
		nw = sf_n(aw);
		
		data = sf_fslice_init(n,nw,sizeof(sf_complex));
		wfld = sf_fslice_init(n,nw,sizeof(sf_complex));
		
		sf_fslice_load(Fd,data,SF_COMPLEX);
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
		sf_setlabel(aw,  "w");
		sf_setunit (aw,"1/s");
		
		amx = sf_iaxa(Fi,1); sf_setlabel(amx,"mx"); sf_oaxa(Fd,amx,1);
		amy = sf_iaxa(Fi,2); sf_setlabel(amy,"my"); sf_oaxa(Fd,amy,2);
		ahx = sf_iaxa(Fi,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fd,ahx,3);
		amz = sf_iaxa(Fi,4); sf_setlabel(amz, "z"); sf_oaxa(Fd,aw ,4);
		
		n = sf_n(amx)*sf_n(amy)*sf_n(ahx);
		
		data = sf_fslice_init(n,nw,sizeof(sf_complex));
		imag = sf_fslice_init(n,nz,sizeof(float));
		
		sf_fslice_load(Fi,imag,SF_FLOAT);
	    } else { /* migration */
		Fd = sf_input ( "in");
		Fi = sf_output("out"); sf_settype(Fi,SF_FLOAT);
		if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
		
		amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fi,amx,1);
		amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fi,amy,2);
		ahx = sf_iaxa(Fd,3); sf_setlabel(ahx,"hx"); sf_oaxa(Fi,ahx,3);
		aw  = sf_iaxa(Fd,4); sf_setlabel(aw , "w"); sf_oaxa(Fi,amz,4);
		
		n  = sf_n(amx)*sf_n(amy)*sf_n(ahx);
		nw = sf_n(aw);

		data = sf_fslice_init(n,nw,sizeof(sf_complex));
		imag = sf_fslice_init(n,nz,sizeof(float));
		
		sf_fslice_load(Fd,data,SF_COMPLEX);
	    } 
	    break;
    }
    
    /*------------------------------------------------------------*/
    cub = camig3_cube(verb,
		      amx,amy,amz,ahx,
		      alx,aly,
		      aw,
		      eps,
		      ompnth,
		      ompchunk);
    
    dsmax = dtmax/cub->amz.d;
    
    /*------------------------------------------------------------*/
    /* init structures */
    tap = taper_init(cub->amx.n,
		     cub->amy.n,
		     cub->ahx.n,
		     SF_MIN(tmx,cub->amx.n-1), 
		     SF_MIN(tmy,cub->amy.n-1), 
		     SF_MIN(thx,cub->ahx.n-1), 
		     true,false,false);

    cam = cam3_init(cub,pmx,pmy,phx,tmx,tmy,thx,dsmax);

    slo = slow3_init(cub,slow,nrmax,dsmax,twoway);
    /*------------------------------------------------------------*/
    weop = camig3_init(cub);
    
    switch(mode[0]) {
	case 'w':
	    cawfl3(weop,cub,cam,tap,slo,inv,data,wfld);
	    break;
	case 'd':
	    cadtm3(weop,cub,cam,tap,slo,inv,data,wfld);
	    break;
	case 'm':
	default:
	    camig3(weop,cub,cam,tap,slo,inv,data,imag);
	    break;
    }

    camig3_close(weop);

    /*------------------------------------------------------------*/
    /* close structures   */
    slow3_close(slo);
    cam3_close(cam);
    taper2d_close(tap);

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    switch(mode[0]) {
	case 'w':
	    sf_fslice_dump(Fw,wfld,SF_COMPLEX);
	    sf_fslice_close(data);
	    sf_fslice_close(wfld);
	    break;
	case 'd':
	    if(inv) sf_fslice_dump(Fd,data,SF_COMPLEX);
	    else    sf_fslice_dump(Fw,wfld,SF_COMPLEX);
	    sf_fslice_close(data);
	    sf_fslice_close(wfld);
	    break;
	case 'm':
	    if(inv) sf_fslice_dump(Fd,data,SF_COMPLEX);
	    else    sf_fslice_dump(Fi,imag,SF_FLOAT);
	    sf_fslice_close(data);
	    sf_fslice_close(imag);
	default:
	    break;
    }
    sf_fslice_close(slow);
    /*------------------------------------------------------------*/


    exit (0);
}
