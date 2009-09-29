/* 3-D zero-offset modeling/migration with extended SSF with time images */

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

#include "zomig3.h"
#include "taper3.h"
#include "ssr3.h"
#include "slow3.h"

#include "weutil.h"
/*^*/

int main (int argc, char *argv[])
{
    char *mode;           /* mode of operation */
    bool verb;            /* verbosity */
    bool inv;             /* forward or adjoint */
    bool causal;          /* causal/anti-causal flag */
    bool twoway;          /* two-way traveltime */
    float eps;            /* dip filter constant */  
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */

    sf_axis amx,amy,amz;
    sf_axis alx,aly;
    sf_axis aw,ae;

    int n,nz,nw;
    float dw,ow;
    
    /* I/O files */
    sf_file Fs=NULL;    /*  slowness file S(nlx,nly,nz   ) */
    sf_file Fi=NULL;    /*     image file R(nmx,nmy,nz   ) */
    sf_file Fd=NULL;    /*      data file D(nmx,nmy,   nw) */
    sf_file Fw=NULL; 

    /* I/O slices */
    fslice slow=NULL;
    fslice imag=NULL;
    fslice data=NULL;
    fslice wfld=NULL;

    int ompchunk=1;
    int ompnth=1;
#ifdef _OPENMP
    int ompath=1; 
#endif

    cub3d cub; /* wavefield hypercube */
    tap3d tap; /* tapering */
    ssr3d ssr; /* SSR operator */
    slo3d slo; /* slowness */

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
    
    /* default mode is migration/modeling */
    if (NULL == (mode = sf_getstring("mode"))) mode = "m";

    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */
    if (!sf_getfloat(  "eps",&eps   ))   eps =  0.01; /* stability parameter */
    if (!sf_getbool(   "inv",&inv   ))   inv = false; /* y=modeling; n=migration */
    if (!sf_getbool("causal",&causal)) causal= false; /* y=causal; n=anti-causal */
    if (!sf_getbool("twoway",&twoway)) twoway=  true; /* two-way traveltime */
    if (!sf_getint(  "nrmax",&nrmax )) nrmax =     1; /* maximum references */
    if (!sf_getfloat("dtmax",&dtmax )) dtmax = 0.004; /* time error */

    if (!sf_getint(    "pmx",&pmx   ))   pmx =     0; /* padding on x */
    if (!sf_getint(    "pmy",&pmy   ))   pmy =     0; /* padding on y*/

    if (!sf_getint(    "tmx",&tmx   ))   tmx =     0; /* taper on x*/
    if (!sf_getint(    "tmy",&tmy   ))   tmy =     0; /* taper on y */

    /* slowness parameters */
    Fs  = sf_input ("slo");
    alx = sf_iaxa(Fs,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs,2); sf_setlabel(aly,"ly");
    amz = sf_iaxa(Fs,3); sf_setlabel(amz,"z" );

    n = sf_n(alx)*sf_n(aly);
    nz = sf_n(amz);

    slow = fslice_init(n,nz,sizeof(float));
    fslice_load(Fs,slow,SF_FLOAT);
    
    switch(mode[0]) {
	case 'w': /* save wavefield */
	    Fd = sf_input ( "in");
	    Fw = sf_output("out"); sf_settype(Fw,SF_COMPLEX);
	    if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
 
	    amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fw,amx,1);
	    amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fw,amy,2);
	    ;                                           sf_oaxa(Fw,amz,3);
	    aw  = sf_iaxa(Fd,3); sf_setlabel(aw ,"w" ); sf_oaxa(Fw,aw ,4);
	    ae  = sf_maxa(1,0,1);

	    n  = sf_n(amx)*sf_n(amy);
	    nw = sf_n(aw);

	    data = fslice_init(n,   nw,sizeof(sf_complex));
	    wfld = fslice_init(n,nz*nw,sizeof(sf_complex));

	    fslice_load(Fd,data,SF_COMPLEX);
	    break;
	case 'd':
	    if (inv) { /*   upward continuation */
		Fw = sf_input ( "in");
		Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
		if (SF_COMPLEX !=sf_gettype(Fw)) sf_error("Need complex data");

		amx = sf_iaxa(Fw,1); sf_setlabel(amx,"mx"); sf_oaxa(Fd,amx,1);
		amy = sf_iaxa(Fw,2); sf_setlabel(amy,"my"); sf_oaxa(Fd,amy,2);
		aw  = sf_iaxa(Fw,3); sf_setlabel(aw , "w"); sf_oaxa(Fd,aw ,3);
		ae  = sf_iaxa(Fw,4); sf_setlabel(ae , "e"); sf_oaxa(Fd,ae ,4);

		n  = sf_n(amx)*sf_n(amy);
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
		aw  = sf_iaxa(Fd,3); sf_setlabel(aw , "w"); sf_oaxa(Fw,aw ,3);
		ae  = sf_iaxa(Fd,4); sf_setlabel(ae , "e"); sf_oaxa(Fw,ae ,4);

		n  = sf_n(amx)*sf_n(amy);
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
		if (!sf_getfloat("ow",&ow)) ow=0.;
		aw = sf_maxa(nw,ow,dw); 
		sf_setlabel(aw,  "w"); 
		sf_setunit (aw,"1/s");
		ae = sf_maxa(1,0,1);

		amx = sf_iaxa(Fi,1); sf_setlabel(amx,"mx"); sf_oaxa(Fd,amx,1);
		amy = sf_iaxa(Fi,2); sf_setlabel(amy,"my"); sf_oaxa(Fd,amy,2);
		amz = sf_iaxa(Fi,3); sf_setlabel(amz, "z"); sf_oaxa(Fd,aw ,3);

		n = sf_n(amx)*sf_n(amy);

		data = fslice_init(n,nw,sizeof(sf_complex));
		imag = fslice_init(n,nz,sizeof(float));

		fslice_load(Fi,imag,SF_FLOAT);			
	    } else { /* migration, output = frequency domain complex cross-correlation at t=0 */
		Fd = sf_input ( "in");
		Fi = sf_output("out"); sf_settype(Fi,SF_COMPLEX);
		if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
		
		amx = sf_iaxa(Fd,1); sf_setlabel(amx,"mx"); sf_oaxa(Fi,amx,1);
		amy = sf_iaxa(Fd,2); sf_setlabel(amy,"my"); sf_oaxa(Fi,amy,2);
		aw  = sf_iaxa(Fd,3); sf_setlabel(aw , "w"); sf_oaxa(Fi,amz,3);
		ae  = sf_maxa(1,0,1); 

		n  = sf_n(amx)*sf_n(amy);
		nw = sf_n(aw);

		data = fslice_init(n,nw,sizeof(sf_complex));
		imag = fslice_init(n,nz,sizeof(sf_complex));
	    
		fslice_load(Fd,data,SF_COMPLEX);
	    }
	    break;
    }
    /*------------------------------------------------------------*/
    cub = zomig3_cube(verb,
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

    slo = slow3_init(cub,slow,nrmax,dsmax,twoway);
    /*------------------------------------------------------------*/
    weop = zomig3_init(cub);

    switch(mode[0]) {
	case 'w':
	    zowfl3(weop,cub,ssr,tap,slo,inv,causal,data,wfld);
	    break;
	case 'd':
	    zodtm3(weop,cub,ssr,tap,slo,inv,causal,data,wfld);
	    break;
	case 'm':
	default:
	    zomig3(weop,cub,ssr,tap,slo,inv,       data,imag);
	    break;
    }

    zomig3_close(weop);

    /*------------------------------------------------------------*/
    /* close structures   */
    slow3_close(slo);
    ssr3_close(ssr);
    taper2d_close(tap);

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
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
	    else    fslice_dump(Fi,imag,SF_COMPLEX);
	    fslice_close(data);
	    fslice_close(imag);
	default:
	    break;
    }
    fslice_close(slow);

    exit (0);
}
