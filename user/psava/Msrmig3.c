/* 3-D S/R migration with extended SSF */

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

#include "srmig3.h"
#include "taper3.h"
#include "ssr3.h"
#include "slow3.h"
#include "img3.h"

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
    bool  cw;             /* converted waves flag */
    char *itype;          /* imaging type 
			     o = zero lag (default)
                             e = extended
			     x = space-lags
			     h = absolute space-lag
			     t = time-lag
			  */
    bool  hsym;
    float vpvs;

    void (*imop)      (cub3d,img3d,int,int);        /* imaging operator apply */
    void (*imop_close)(cub3d,img3d,sf_fslice,sf_fslice);  /* imaging operator close */

    sf_axis amx,amy,amz;
    sf_axis alx,aly;
    sf_axis aw,ae;
    sf_axis ahx,ahy,ahz,aht;
    sf_axis acx,acy,acz;
    int     jcx,jcy,jcz;
    sf_axis ahh,aha,ahb; /* |h|,alpha,beta */

    int             nhx,nhy,nhz,nw;
    int n,nz,nx,ny, nht,nh,nha,nhb;
    float           oht,oh,oha,ohb;
    float dz,dx,dy, dht,dh,dha,dhb;

    /* I/O files */
    sf_file Fw_s=NULL,Fw_r=NULL; /* wavefield file W ( nx, ny,nw) */
    sf_file Fs_s=NULL,Fs_r=NULL; /*  slowness file S (nlx,nly,nz) */
    sf_file Fi=NULL;             /*     image file R ( nx, ny,nz) */
    sf_file Fc=NULL;             /*      cigs file C */

    /* I/O slices */
    sf_fslice wfl_s=NULL,wfl_r=NULL;
    sf_fslice slo_s=NULL,slo_r=NULL;
    sf_fslice imag=NULL;
    sf_fslice cigs=NULL;

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
    img3d img; /* imaging  */

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

    if (NULL == (itype = sf_getstring("itype"))) itype = "o";
    /* imaging condition type
       o = zero lag (default)
       e = extended
       x = space-lags
       h = space-lags magnitude
       t = time-lag
    */

    if (!sf_getbool(  "verb",&verb ))  verb =  true; /* verbosity flag */
    if (!sf_getfloat(  "eps",&eps  ))   eps =  0.01; /* stability parameter */
    if (!sf_getbool("twoway",&twoway)) twoway=false; /* two-way traveltime */
    if (!sf_getint(  "nrmax",&nrmax)) nrmax =     1; /* max number of refs */
    if (!sf_getfloat("dtmax",&dtmax)) dtmax = 0.004; /* max time error */

    if (!sf_getint(    "pmx",&pmx  ))   pmx =     0; /* padding on x */
    if (!sf_getint(    "pmy",&pmy  ))   pmy =     0; /* padding on y */

    if (!sf_getint(    "tmx",&tmx  ))   tmx =     0; /* taper on x   */
    if (!sf_getint(    "tmy",&tmy  ))   tmy =     0; /* taper on y   */

    if (!sf_getfloat( "vpvs",&vpvs))   vpvs =    1.; /* Vp/Vs ratio */

    /*------------------------------------------------------------*/
    /* slowness parameters */
    ;      Fs_s = sf_input("slo");
    if(cw) Fs_r = sf_input("sls");
    else   Fs_r = sf_input("slo");

    alx = sf_iaxa(Fs_s,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Fs_s,2); sf_setlabel(aly,"ly");
    amz = sf_iaxa(Fs_s,3); sf_setlabel(amz,"z" );
    /* test here if slo and sls have similar sizes */

    n = sf_n(alx)*sf_n(aly);
    nz = sf_n(amz);
    dz = sf_d(amz);

    slo_s = sf_fslice_init(n, nz, sizeof(float));
    slo_r = sf_fslice_init(n, nz, sizeof(float));
    sf_fslice_load(Fs_s,slo_s,SF_FLOAT);
    sf_fslice_load(Fs_r,slo_r,SF_FLOAT);

    /*------------------------------------------------------------*/    
    /* WAVEFIELD/IMAGE */

    Fw_s = sf_input( "in"); /*   source data */
    Fw_r = sf_input("rwf"); /* receiver data */   
    if (SF_COMPLEX != sf_gettype(Fw_s)) sf_error("need complex   source data");
    if (SF_COMPLEX != sf_gettype(Fw_r)) sf_error("need complex receiver data");
    
    amx = sf_iaxa(Fw_s,1); sf_setlabel(amx,"x");
    amy = sf_iaxa(Fw_s,2); sf_setlabel(amy,"y");
    aw  = sf_iaxa(Fw_s,3); sf_setlabel(aw, "w");
    ae  = sf_iaxa(Fw_s,4); sf_setlabel(ae, "e"); /* experiments */

    nx=sf_n(amx); dx=sf_d(amx);
    ny=sf_n(amy); dy=sf_d(amy);

    Fi   = sf_output("out"); sf_settype(Fi,SF_FLOAT);
    sf_oaxa(Fi,amx,1);
    sf_oaxa(Fi,amy,2);
    sf_oaxa(Fi,amz,3);
    sf_putint(Fi,"n4",1);
    sf_putint(Fi,"n5",1);

    imag = sf_fslice_init( nx*ny*nz,1,sizeof(float));

    /*------------------------------------------------------------*/
    /* CIGS */

    if(!sf_getint ("jcx",&jcx) || nx==1) jcx=1;
    if(!sf_getint ("jcy",&jcy) || ny==1) jcy=1;
    if(!sf_getint ("jcz",&jcz) || nz==1) jcz=1;
    /* CIGs windowing */

    acx = sf_maxa(SF_MAX(1,nx/jcx),sf_o(amx),dx*jcx); sf_setlabel(acx,"cx");
    acy = sf_maxa(SF_MAX(1,ny/jcy),sf_o(amy),dy*jcy); sf_setlabel(acy,"cy");
    acz = sf_maxa(SF_MAX(1,nz/jcz),sf_o(amz),dz*jcz); sf_setlabel(acz,"cz");
    n = sf_n(acx)*sf_n(acy)*sf_n(acz);

    Fc = sf_output("cig");
    /* Output file with Common Image Gathers */
    sf_settype(Fc,SF_FLOAT);
    sf_oaxa(Fc,acx,1);
    sf_oaxa(Fc,acy,2);
    sf_oaxa(Fc,acz,3);

    /* A file with CIGs will be written even when itype=o and CIGs are the same
    as the image. This is in order to allow switching between imaging
    conditions without changing SConstruct files. Changing this will break
    several SConstruct files in RSFSRC/book . */

    /*------------------------------------------------------------*/
    /* wavefield hypercube */
    cub = srmig3_cube(verb,
		      amx,amy,amz,
		      alx,aly,
		      aw,
		      ae,
		      eps,
		      ompnth,
		      ompchunk);
    
    dsmax = dtmax/cub->amz.d;

    /*------------------------------------------------------------*/
    /* init output files */
    switch(itype[0]) {
	case 'e':
	    if(verb) sf_warning("E.I.C.");

	    /* x-lags */
	    if(!sf_getint("nhx",&nhx) || nx==1) nhx=1;
	    if(!sf_getint("nhy",&nhy) || ny==1) nhy=1;
	    if(!sf_getint("nhz",&nhz) || nz==1) nhz=1;
	    ahx = sf_maxa(nhx,0.,dx); sf_setlabel(ahx,"hx");
	    ahy = sf_maxa(nhy,0.,dy); sf_setlabel(ahy,"hy");
	    ahz = sf_maxa(nhz,0.,dz); sf_setlabel(ahz,"hz");
	    if(!sf_getbool("hsym",&hsym)) hsym = false;
	    if(hsym) {
		if(nhx>1) { sf_seto(ahx,-nhx*dx); sf_setn(ahx,nhx*2); }
		if(nhy>1) { sf_seto(ahy,-nhy*dy); sf_setn(ahy,nhy*2); }
		if(nhz>1) { sf_seto(ahz,-nhz*dz); sf_setn(ahz,nhz*2); }
	    }

	    /* t-lag */
	    if(!sf_getint  ("nht",&nht)) nht=1;
	    if(!sf_getfloat("oht",&oht)) oht=0;
	    if(!sf_getfloat("dht",&dht)) dht=0.1;
	    aht = sf_maxa(nht,oht,dht); sf_setlabel(aht,"ht");

	    sf_oaxa(Fc,ahx,4); sf_raxa(ahx);
	    sf_oaxa(Fc,ahy,5); sf_raxa(ahy);
	    sf_oaxa(Fc,ahz,6); sf_raxa(ahz);
	    sf_oaxa(Fc,aht,7); sf_raxa(aht);

	    sf_raxa(acx);
	    sf_raxa(acy);
	    sf_raxa(acz);

	    cigs = sf_fslice_init(n,
			       sf_n(ahx)*sf_n(ahy)*sf_n(ahz)*sf_n(aht),
			       sizeof(float));

	    img=img3e_init(cub,imag,cigs,jcx,jcy,jcz,ahx,ahy,ahz,aht);
	    imop       = img3e;
	    imop_close = img3e_close;
	    break;
	case 't':
	    if(verb) sf_warning("t-lag I.C.");
	    if(!sf_getint  ("nht",&nht)) nht=1;
	    if(!sf_getfloat("oht",&oht)) oht=0;
	    if(!sf_getfloat("dht",&dht)) dht=0.1;
	    aht = sf_maxa(nht,oht,dht); sf_setlabel(aht,"ht");
	    sf_oaxa(Fc,aht,4);

	    cigs = sf_fslice_init(n,
			       sf_n(aht),
			       sizeof(float));

	    img=img3t_init(cub,imag,cigs,jcx,jcy,jcz,aht);
	    imop       = img3t;
	    imop_close = img3t_close;
	    break;
	case 'x':
	    if(verb) sf_warning("x-lags I.C.");
	    if(!sf_getint("nhx",&nhx) || nx==1) nhx=1;
	    if(!sf_getint("nhy",&nhy) || ny==1) nhy=1;
	    if(!sf_getint("nhz",&nhz) || nz==1) nhz=1;
	    ahx = sf_maxa(nhx,0.,dx); sf_setlabel(ahx,"hx");
	    ahy = sf_maxa(nhy,0.,dy); sf_setlabel(ahy,"hy");
	    ahz = sf_maxa(nhz,0.,dz); sf_setlabel(ahz,"hz");

	    if(!sf_getbool("hsym",&hsym)) hsym = false;
	    if(hsym) {
		if(nhx>1) { sf_seto(ahx,-nhx*dx); sf_setn(ahx,nhx*2); }
		if(nhy>1) { sf_seto(ahy,-nhy*dy); sf_setn(ahy,nhy*2); }
		if(nhz>1) { sf_seto(ahz,-nhz*dz); sf_setn(ahz,nhz*2); }
	    }
	    sf_oaxa(Fc,ahx,4); sf_raxa(ahx);
	    sf_oaxa(Fc,ahy,5); sf_raxa(ahy);
	    sf_oaxa(Fc,ahz,6); sf_raxa(ahz);

	    cigs = sf_fslice_init(n,
			       sf_n(ahx)*sf_n(ahy)*sf_n(ahz),
			       sizeof(float));

	    img=img3x_init(cub,imag,cigs,jcx,jcy,jcz,ahx,ahy,ahz);
	    imop       = img3x;
	    imop_close = img3x_close;
	    break;
	case 'h':
	    if(verb) sf_warning("|x|-lag I.C.");

	    if(!sf_getint  ("nhh",&nh)) nh=1;
	    if(!sf_getfloat("ohh",&oh)) oh=0;
	    if(!sf_getfloat("dhh",&dh)) dh=0.1;
	    ahh = sf_maxa(nh,oh,dh); sf_setlabel(ahh,"hh");
	    
	    /* longitude */
	    if(!sf_getint  ("nha",&nha)) nha=180;
	    if(!sf_getfloat("oha",&oha)) oha=0;   
	    if(!sf_getfloat("dha",&dha)) dha=2.0; 
        oha *= SF_PI/180;
        dha *= SF_PI/180;
	    if(nz==1) { nha=1; oha=0.;       dha=SF_PI;}
	    aha = sf_maxa(nha,oha,dha); sf_setlabel(aha,"ha");

	    /* latitude */
	    if(!sf_getint  ("nhb",&nhb)) nhb=180;
	    if(!sf_getfloat("ohb",&ohb)) ohb=0;   
	    if(!sf_getfloat("dhb",&dhb)) dhb=2.0; 
        ohb *= SF_PI/180;
        dhb *= SF_PI/180;
	    if(nx==1) { nhb=1; ohb=SF_PI/2.; dhb=SF_PI;}
	    if(ny==1) { nhb=1; ohb=0.;       dhb=SF_PI;}
	    ahb = sf_maxa(nhb,ohb,dhb); sf_setlabel(ahb,"hb");

	    sf_raxa(ahh);
	    sf_raxa(aha);
	    sf_raxa(ahb);

	    sf_oaxa(Fc,ahh,4);
	    cigs = sf_fslice_init(n,
			       nh,
			       sizeof(float));

	    img=img3h_init(cub,imag,cigs,jcx,jcy,jcz,ahh,aha,ahb,vpvs);
	    imop       = img3h;
	    imop_close = img3h_close;
	    break;
	case 'o':
	default:
	    if(verb) sf_warning("C.I.C.");
	    cigs = sf_fslice_init(n,1,sizeof(float));
	    ahx = sf_maxa(1,0.,0); sf_setlabel(ahx,"");
	    sf_oaxa(Fc,ahx,4);

	    img=img3o_init(cub,imag,cigs,jcx,jcy,jcz);
	    imop       = img3o;
	    imop_close = img3o_close;
	    break;
    }

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    nw = sf_n(aw)*sf_n(ae);

    wfl_s = sf_fslice_init( nx*ny, nw,sizeof(sf_complex));
    wfl_r = sf_fslice_init( nx*ny, nw,sizeof(sf_complex));

    sf_fslice_load(Fw_s,wfl_s,SF_COMPLEX);
    sf_fslice_load(Fw_r,wfl_r,SF_COMPLEX);

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
    /* MIGRATION */
    weop = srmig3_init(cub);

    srmig3(weop,  /* shot-record migration operator */
	   cub,   /* wavefield hypercube dimensions */
	   ssr,   /* SSR operator */
	   tap,   /* tapering operator */
	   s_s,   /* source slowness */
	   s_r,   /* receiver slowness */
	   img,   /* image */
	   wfl_s, /* source wavefield */
	   wfl_r, /* receiver wavefield */
	   imag,  /* image */
	   cigs,  /* CIGs */
	   imop   /* imaging operator */
	);

    srmig3_close(weop);

    /*------------------------------------------------------------*/
    /* close structures */
    slow3_close(s_s);
    slow3_close(s_r);
    ssr3_close(ssr);
    taper2d_close(tap);
    
    if(verb) sf_warning("imop close");
    imop_close(cub,img,imag,cigs); 

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    if(verb) sf_warning("dump imag");
    sf_fslice_dump(Fi,imag,SF_FLOAT);

    if(verb) sf_warning("dump cigs");
    sf_fslice_dump(Fc,cigs,SF_FLOAT);

    /*------------------------------------------------------------*/
    sf_fslice_close(slo_s);
    sf_fslice_close(slo_r);
    sf_fslice_close(wfl_s);
    sf_fslice_close(wfl_r);
    sf_fslice_close(imag);
    sf_fslice_close(cigs);

    /*------------------------------------------------------------*/


    exit (0);
}
