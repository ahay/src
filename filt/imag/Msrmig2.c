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
#include "srmig2.h"
#include "img2.h"

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
			     o = zero lag (default)
			     x = space-lags
			     h = absolute space-lag
			     t = time-lag
			  */
    bool  hsym;
    float vpvs;
/*    int dec;*/
/*    float deceps;*/
    void (*imop)(int);
    void (*close)(fslice,fslice);

    sf_axis amx,amy,amz;
    sf_axis alx,aly;
    sf_axis aw,ae;
    sf_axis ahx,ahy,ahz,aht;
    sf_axis acx,acy,acz;
    int     jcx,jcy,jcz;
    sf_axis ahh,aha,ahb; /* |h|,alpha,beta */

    int n,nz,nx,ny,nh,nhx,nhy,nhz,nha,nhb,nw;
    float dz,dx,dy,oh,dh,oha,dha,ohb,dhb;

    /* I/O files */
    sf_file Fw_s=NULL,Fw_r=NULL; /* wavefield file W ( nx, ny,nw) */
    sf_file Fs_s=NULL,Fs_r=NULL; /*  slowness file S (nlx,nly,nz) */
    sf_file Fi=NULL;             /*     image file R ( nx, ny,nz) */
    sf_file Fc=NULL;             /*      cigs file C */

    /* I/O slices */
    fslice wfl_s=NULL,wfl_r=NULL;
    fslice slo_s=NULL,slo_r=NULL;
    fslice imag=NULL;
    fslice cigs=NULL;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (NULL == (itype = sf_getstring("itype"))) itype = "o";
    /* imaging condition type 
       o = zero lag (default)
       x = space-lags
       h = space-lags magnitude
       t = time-lag
    */

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

    if (!sf_getfloat( "vpvs",&vpvs))   vpvs = 1.;
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
    dz = sf_d(amz);

    ;      slo_s = fslice_init(n, nz, sizeof(float));
    if(cw) slo_r = fslice_init(n, nz, sizeof(float));
    ;      fslice_load(Fs_s,slo_s,SF_FLOAT);
    if(cw) fslice_load(Fs_r,slo_r,SF_FLOAT);

    /*------------------------------------------------------------*/    
    /* WAVEFIELD/IMAGE */

    Fw_s = sf_input( "in");
    Fw_r = sf_input("rwf");
    
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

    imag = fslice_init( nx*ny*nz,1,sizeof(float));
    
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

    Fc = sf_output("cig"); sf_settype(Fc,SF_FLOAT);
    sf_oaxa(Fc,acx,1);
    sf_oaxa(Fc,acy,2);
    sf_oaxa(Fc,acz,3);

    /*------------------------------------------------------------*/
    /* init output files */

/*    if (!sf_getint("dec",&dec)) dec=0;*/
/*    if (dec) {*/
/*	if (!sf_getfloat("deceps",&deceps)) deceps=0.01;*/
/*	img2dec(dec,deceps);*/
/*    }*/

    switch(itype[0]) {
	case 't':
	    if(verb) sf_warning("time-lag I.C.");
	    if(!sf_getint  ("nht",&nh)) nh=1;
	    if(!sf_getfloat("oht",&oh)) oh=0;
	    if(!sf_getfloat("dht",&dh)) dh=0.1;
	    aht = sf_maxa(nh,oh,dh); sf_setlabel(aht,"ht");
	    
	    sf_oaxa(Fc,aht,4);
	    cigs = fslice_init(n*nh,1,sizeof(float));

	    img2t_init(amx,amy,amz,jcx,jcy,jcz,aht,aw,imag);
	    imop  = img2t;
	    close = img2t_close;
	    break;
	case 'x':
	    if(verb) sf_warning("space-lags I.C.");
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

	    n *= sf_n(ahx)*sf_n(ahy)*sf_n(ahz);
	    cigs = fslice_init(n,1,sizeof(float));

	    img2x_init(amx,amy,amz,jcx,jcy,jcz,ahx,ahy,ahz,imag);
	    imop  = img2x;
	    close = img2x_close;
	    break;
	case 'g':
	case 'h':
	    if(verb) sf_warning("space-lag magnitude I.C.");

	    if(!sf_getint  ("nhh",&nh)) nh=1;
	    if(!sf_getfloat("ohh",&oh)) oh=0;
	    if(!sf_getfloat("dhh",&dh)) dh=0.1;
	    ahh = sf_maxa(nh,oh,dh); sf_setlabel(ahh,"hh");
	    
	    /* longitude */
	    if(!sf_getint  ("nha",&nha)) nha=180;
	    if(!sf_getfloat("oha",&oha)) oha=0;   oha *= SF_PI/180;
	    if(!sf_getfloat("dha",&dha)) dha=2.0; dha *= SF_PI/180;
	    if(nz==1) { nha=1; oha=0.;       dha=SF_PI;}
	    aha = sf_maxa(nha,oha,dha); sf_setlabel(aha,"ha");

	    /* latitude */
	    if(!sf_getint  ("nhb",&nhb)) nhb=180;
	    if(!sf_getfloat("ohb",&ohb)) ohb=0;   ohb *= SF_PI/180;
	    if(!sf_getfloat("dhb",&dhb)) dhb=2.0; dhb *= SF_PI/180;
	    if(nx==1) { nhb=1; ohb=SF_PI/2.; dhb=SF_PI;}
	    if(ny==1) { nhb=1; ohb=0.;       dhb=SF_PI;}
	    ahb = sf_maxa(nhb,ohb,dhb); sf_setlabel(ahb,"hb");

	    sf_raxa(ahh);
	    sf_raxa(aha);
	    sf_raxa(ahb);

	    sf_oaxa(Fc,ahh,4);
	    cigs = fslice_init( n*nh,1,sizeof(float));

	    img2h_init(amx,amy,amz,jcx,jcy,jcz,ahh,aha,ahb,aw,imag,vpvs);
	    imop  = img2h;
	    close = img2h_close;
	    break;
	case 'o':
	default:
	    if(verb) sf_warning("zero-lag I.C.");
	    cigs = fslice_init(n,1,sizeof(float));

	    img2o_init(amx,amy,amz,jcx,jcy,jcz,imag);
	    imop  = img2o;
	    close = img2o_close;
	    break;
    }

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    nw = sf_n(aw)*sf_n(ae);

    wfl_s = fslice_init( nx*ny, nw,sizeof(sf_complex));
    wfl_r = fslice_init( nx*ny, nw,sizeof(sf_complex));

    fslice_load(Fw_s,wfl_s,SF_COMPLEX);
    fslice_load(Fw_r,wfl_r,SF_COMPLEX);

    /*------------------------------------------------------------*/
    /* MIGRATION */
    srmig2_init (verb,eps,dtmax,
		 ae,aw,amx,amy,amz,alx,aly,
		 tmx,tmy,pmx,pmy);
    
    if(cw) { 
	srmig2_cw_init (dtmax,nrmax,slo_s,slo_r);
	srmig2_cw(wfl_s,wfl_r,imag,cigs, imop); 
	srmig2_cw_close();
    } else { 
	srmig2_pw_init (dtmax,nrmax,slo_s);
	srmig2_pw(wfl_s,wfl_r,imag,cigs, imop); 
	srmig2_pw_close();
    }
    
    srmig2_close();
    close(imag,cigs); 

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    if(verb) sf_warning("dump imag");
    fslice_dump(Fi,imag,SF_FLOAT);

    if(verb) sf_warning("dump cigs");
    fslice_dump(Fc,cigs,SF_FLOAT)

    ;      fslice_close(slo_s);
    if(cw) fslice_close(slo_r);
    ;      fslice_close(wfl_s);
    ;      fslice_close(wfl_r);
    ;      fslice_close(imag);

    exit (0);
}
