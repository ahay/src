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
			     o = zero offset (default)
			     x = space offset
			     t = time offset
			  */
    bool hsym;

    axa amx,amy,amz;
    axa alx,aly;
    axa aw,ae,aj;
    axa ahx,ahy,ahz,aht;
    axa acx,acy,acz;
    int jcx,jcy,jcz;

    /* I/O files */
    sf_file Fw_s,Fw_r; /* wavefield file W ( nx, ny,nw) */
    sf_file Fs_s,Fs_r; /*  slowness file S (nlx,nly,nz) */
    sf_file Fi;        /*     image file R ( nx, ny,nz) */
    sf_file Fc;        /*      cigs file C */

    /* I/O slices */
    fslice wfl_s,wfl_r;
    fslice slo_s,slo_r;
    fslice imag;
    fslice cigs;

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
    
    if (SF_COMPLEX != sf_gettype(Fw_s)) sf_error("Need complex   source data");
    if (SF_COMPLEX != sf_gettype(Fw_r)) sf_error("Need complex receiver data");

    aj.n=1; aj.o=0; aj.d=1; aj.l=" ";
    
    iaxa(Fw_s,&amx,1); amx.l="mx";
    iaxa(Fw_s,&amy,2); amy.l="my";
    iaxa(Fw_s,&aw,3);  aw.l = "w";
    iaxa(Fw_s,&ae,4);  ae.l = "e"; /* experiments */

    Fi   = sf_output("out"); sf_settype(Fi,SF_FLOAT);
    oaxa(Fi,&amx,1);
    oaxa(Fi,&amy,2);
    oaxa(Fi,&amz,3);
    oaxa(Fi,&aj, 4);
    oaxa(Fi,&aj, 5);

    imag = fslice_init( amx.n*amy.n*amz.n,1,sizeof(float));    
    /*------------------------------------------------------------*/
    /* CIGS */

    if(!sf_getint  ("jcx",&jcx)) jcx=1;
    if(!sf_getint  ("jcy",&jcy)) jcy=1;
    if(!sf_getint  ("jcz",&jcz)) jcz=1;

    acx.o = amx.o;
    acx.n = amx.n / jcx;
    acx.d = amx.d * jcx;
    acx.l="cx";
    
    acy.o = amy.o;
    acy.n = amy.n / jcy;
    acy.d = amy.d * jcy;
    acy.l="cy";
    
    acz.o = amz.o;
    acz.n = amz.n / jcz;
    acz.d = amz.d * jcz;
    acz.l="cz";	    

    Fc = sf_output("cig"); sf_settype(Fc,SF_FLOAT);
    oaxa(Fc,&acx,1);
    oaxa(Fc,&acy,2);
    oaxa(Fc,&acz,3);

    /*------------------------------------------------------------*/
    /* init output files */

    switch(itype[0]) {
	case 't': /* time offset imaging condition */
	    if(!sf_getint  ("nht",&aht.n)) aht.n=1;
	    if(!sf_getfloat("oht",&aht.o)) aht.o=0;
	    if(!sf_getfloat("dht",&aht.d)) aht.d=0.1;
	    aht.l="ht";
	    
	    oaxa(Fc,&aht,4);
	    cigs = fslice_init( acx.n*acy.n*acz.n*aht.n,1,sizeof(float));

	    img2t_init(amx,amy,amz,jcx,jcy,jcz,aht,aw,imag);
	    break;
	case 'x': /* space offset imaging condition */
	    if(!sf_getint("nhx",&ahx.n)) ahx.n=1;
	    if(!sf_getint("nhy",&ahy.n)) ahy.n=1;
	    if(!sf_getint("nhz",&ahz.n)) ahz.n=1;
	    ahx.o=0;     ahy.o=0;     ahz.o=0.;
	    ahx.d=amx.d; ahy.d=amy.d; ahz.d=amz.d;
	    ahx.l="hx";  ahy.l="hy";  ahz.l="hz";

	    if(!sf_getbool("hsym",&hsym)) hsym = false;
	    if(hsym) {
		if(ahx.n>1) { ahx.o = - ahx.n * ahx.d; ahx.n *=2; }
		if(ahy.n>1) { ahy.o = - ahy.n * ahy.d; ahy.n *=2; }
		if(ahz.n>1) { ahz.o = - ahz.n * ahz.d; ahz.n *=2; }
	    }

	    oaxa(Fc,&ahx,4);
	    oaxa(Fc,&ahy,5);
	    oaxa(Fc,&ahz,6);
	    cigs = fslice_init( acx.n*acy.n*acz.n*ahx.n*ahy.n*ahz.n,1,sizeof(float));

	    img2x_init(amx,amy,amz,jcx,jcy,jcz,ahx,ahy,ahz,imag);
	    break;
	case 'o': /* zero offset imaging condition */
	default:

	    cigs = fslice_init( acx.n*acy.n*acz.n,1,sizeof(float));

	    img2o_init(amx,amy,amz,jcx,jcy,jcz,imag);
	    break;
    }

    /* slice management (temp files) */
    wfl_s = fslice_init( amx.n*amy.n, aw.n*ae.n,sizeof(float complex));
    wfl_r = fslice_init( amx.n*amy.n, aw.n*ae.n,sizeof(float complex));

    fslice_load(Fw_s,wfl_s,SF_COMPLEX);
    fslice_load(Fw_r,wfl_r,SF_COMPLEX);

    /*------------------------------------------------------------*/
    /* MIGRATION */

    srmig2_init (verb,eps,dtmax,
		 ae,aw,amx,amy,amz,alx,aly,
		 tmx,tmy,pmx,pmy);
    
    if(cw) { 
	srmig2_cw_init (dtmax,nrmax,slo_s,slo_r);
	switch(itype[0]) {
	    case 't':          srmig2_cw(wfl_s,wfl_r,imag,cigs, &img2t); break;
	    case 'x':          srmig2_cw(wfl_s,wfl_r,imag,cigs, &img2x); break;
	    case 'o': default: srmig2_cw(wfl_s,wfl_r,imag,cigs, &img2o); break;
	}
	srmig2_cw_close();
    } else { 
	srmig2_pw_init (dtmax,nrmax,slo_s);
	switch(itype[0]) {
	    case 't':          srmig2_pw(wfl_s,wfl_r,imag,cigs, &img2t); break;
	    case 'x':          srmig2_pw(wfl_s,wfl_r,imag,cigs, &img2x); break;
	    case 'o': default: srmig2_pw(wfl_s,wfl_r,imag,cigs, &img2o); break;
	}
	srmig2_pw_close();
    }
    
    srmig2_close();

    switch(itype[0]) {
	case 't':          img2t_close(imag,cigs); break;
	case 'x':          img2x_close(imag,cigs); break;
	case 'o': default: img2o_close(imag,cigs); break;
    }

    /*------------------------------------------------------------*/
    /* slice management (temp files) */
    fslice_dump(Fi,imag,SF_FLOAT);
    fslice_dump(Fc,cigs,SF_FLOAT)

    ;      fslice_close(slo_s);
    if(cw) fslice_close(slo_r);
    ;      fslice_close(wfl_s);
    ;      fslice_close(wfl_r);
    ;      fslice_close(imag);

    exit (0);
}
