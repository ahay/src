/* 3-D modeling/migration with extended SSF */

/*
  Copyright (C) 2010 Colorado School of Mines

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

#include "wei.h"
#include "weitap.h"
#include "weissr.h"
#include "weislo.h"
#include "weiutl.h"
/*^*/

int main (int argc, char *argv[])
{
    bool  verb;         /* verbosity */
    int ompnth=1;       /* number of threads */
    char *irun;         /* I run... */

    sf_file Fslo=NULL;  /*  slowness file  S(nlx,nly,nz   ) */
    sf_file Fsou=NULL;  /*      data file Ds(nmx,nmy,   nw) */
    sf_file Frec=NULL;  /*      data file Dr(nmx,nmy,   nw) */
    sf_file Fwfl=NULL;  /* wavefield file  W(nmx,nmy,nz,nw) */
    sf_file Fcic=NULL;  /* CIC image file  R(nmx,nmy,nz   ) */
    sf_file Feic=NULL;  /* EIC image file  R(hx,hy,hz,tt,c) */
    sf_file Fcoo=NULL;  /*     coord file  C(nc) */

    weicub3d cub;       /* hypercube */
    weitap3d tap;       /* tapering */
    weissr3d ssr;       /* SSR operator */
    weislo3d slo;       /* slowness */
    weiop3d weop;       /* WEI operator */
    weiop3f weof;       /* WEI operator */
    weico3d eico;       /* eic coordinates */

    bool  causal;       /* causality flag */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

#ifdef _OPENMP
    ompnth=omp_init(); /* OMP parameters */
#endif

    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */
    if (NULL == (irun = sf_getstring("irun"))) irun = "w";

    if(verb) fprintf(stderr,"init cube...");
    cub = wei_cube(verb,ompnth);
    if(verb) fprintf(stderr,"OK\n");

    /*------------------------------------------------------------*/
    Fslo = sf_input ("slo");
    if (SF_FLOAT !=sf_gettype(Fslo)) sf_error("Need float slowness");
    weislo_inp(cub,Fslo); /* input slowness */

    /*------------------------------------------------------------*/
    switch(irun[0]) {

	case 'e': /* EIC  hx-hy-hz */
	case 'h': /* HIC: hx-hy-1  */
	    if(verb) sf_warning("EIC OUT");

	    Fsou = sf_input ( "in");
	    if (SF_COMPLEX !=sf_gettype(Fsou)) sf_error("Need complex sdat");
	    Frec = sf_input ("dat");
	    if (SF_COMPLEX !=sf_gettype(Frec)) sf_error("Need complex rdat");

	    weidat_inp(cub,Fsou); /* input data dimensions */

	    /* CIC */
	    Fcic = sf_output("out"); sf_settype(Fcic,SF_FLOAT);
	    weicic_out(cub,Fcic); /* prepare output image */

	    /* EIC */
	    Fcoo = sf_input ("coo");
	    if (SF_FLOAT !=sf_gettype(Fcoo)) sf_error("Need float coordinates");

	    if(irun[0]=='e')
		weieic_inp(cub,Fcoo); /* input coordinates */
	    else
		weihic_inp(cub,Fcoo); /* input coordinates */

	    Feic = sf_output("cip"); sf_settype(Feic,SF_FLOAT);
	    weieic_out(cub,Feic); /* prepare output image */

	    break;

	case 'c': /* CIC */
	    if(verb) sf_warning("CIC OUT");

	    Fsou = sf_input ( "in");
	    if (SF_COMPLEX !=sf_gettype(Fsou)) sf_error("Need complex sdat");
	    Frec = sf_input ("dat");
	    if (SF_COMPLEX !=sf_gettype(Frec)) sf_error("Need complex rdat");

	    weidat_inp(cub,Fsou); /* input data dimensions */

	    Fcic = sf_output("out"); sf_settype(Fcic,SF_FLOAT);
	    weicic_out(cub,Fcic); /* prepare output image */
	    break;
	case 'f': /* CIC */
	    if(verb) sf_warning("FCIC OUT");

	    Fsou = sf_input ( "in");
	    if (SF_COMPLEX !=sf_gettype(Fsou)) sf_error("Need complex sdat");
	    Frec = sf_input ("dat");
	    if (SF_COMPLEX !=sf_gettype(Frec)) sf_error("Need complex rdat");

	    weidat_inp(cub,Fsou); /* input data dimensions */

	    Fcic = sf_output("out"); sf_settype(Fcic,SF_FLOAT);
	    weific_out(cub,Fcic); /* prepare output image */
	    break;
	case 'd': /* DTM */
	    if(verb) sf_warning("DTM OUT");

	    Fsou = sf_input ( "in");
	    if (SF_COMPLEX !=sf_gettype(Fsou)) sf_error("Need complex data");

	    weidat_inp(cub,Fsou); /* input data dimensions */

	    Fwfl = sf_output("out"); sf_settype(Fwfl,SF_COMPLEX);
	    weidat_out(cub,Fwfl); /* prepare output datumed data */

	    break;

        case 'w': /* wavefield */
	default:
	    if(verb) sf_warning("WFL OUT");

	    Fsou = sf_input ( "in");
	    if (SF_COMPLEX !=sf_gettype(Fsou)) sf_error("Need complex data");

	    weidat_inp(cub,Fsou); /* input data dimensions */

	    Fwfl = sf_output("out"); sf_settype(Fwfl,SF_COMPLEX);
	    weiwfl_out(cub,Fwfl); /* prepare output wavefield */
	    break;
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(verb) wei_report(cub);
    /*------------------------------------------------------------*/
       
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"init slo...");
    slo = weislo_init(cub,Fslo);
    if(verb) fprintf(stderr,"OK\n");
    if(verb) weislo_report(slo);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"init tap...");
    tap = weitap_init(cub);
    if(verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"init ssr...");
    ssr = weissr_init(cub,slo);
    if(verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    switch(irun[0]) {

	case 'e': /* extended IC: hx-hy-hz */

	    weop = weieic_init(cub);
	    eico = weicoo_init(cub,Fcoo);

	    weieic(weop,cub,ssr,tap,slo,eico,Fsou,Frec,Fcic,Feic);

	    weicoo_close(eico);
	    weieic_close(weop);
	    
	    break;
	case 'h': /* extended IC: hx-hy-1 */

	    weop = weihic_init(cub);
	    eico = weicoo_init(cub,Fcoo);

	    weihic(weop,cub,ssr,tap,slo,eico,Fsou,Frec,Fcic,Feic);

	    weicoo_close(eico);
	    weihic_close(weop);
	    
	    break;
	case 'c': /* conventional IC */

	    weop = weicic_init(cub);
	    weicic(weop,cub,ssr,tap,slo,Fsou,Frec,Fcic);
	    weicic_close(weop);

	    break;
	case 'f': /* conventional IC */

	    weof = weific_init(cub);
	    weific(weof,cub,ssr,tap,slo,Fsou,Frec,Fcic);
	    weific_close(weof);
	    break;
	case 'd': /* datuming */

	    if (!sf_getbool("causal",&causal)) causal=false; /* causality */
	    weop = weidtm_init(cub);
	    weidtm(weop,cub,ssr,tap,slo,Fsou,Fwfl,causal);
	    weidtm_close(weop);

	    break;
	case 'w': /* wavefield */
	default:

	    if (!sf_getbool("causal",&causal)) causal=false; /* causality */
	    weop = weiwfl_init(cub);
	    weiwfl(weop,cub,ssr,tap,slo,Fsou,Fwfl,causal);
	    weiwfl_close(weop);
	    break;
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"close structures...");
    weislo_close(slo);
    weissr_close(cub,ssr);
    weitap_close(tap);
    if(verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    exit (0);
}
