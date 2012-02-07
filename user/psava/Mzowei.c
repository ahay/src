/* 3-D zero-offset modeling/migration */

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
    bool   adj;         /* adjoint flag */
    int ompnth=1;       /* number of threads */

    sf_file Fslo=NULL;  /*  slowness file  S(nlx,nly,nz   ) */
    sf_file Fdat=NULL;  /*      data file Dr(nmx,nmy,   nw) */
    sf_file Fsou=NULL;  /*      data file Dr(nmx,nmy,   nw) */
    sf_file Fcic=NULL;  /* CIC image file  R(nmx,nmy,nz   ) */

    weicub3d cub;       /* hypercube */
    weitap3d tap;       /* tapering */
    weissr3d ssr;       /* SSR operator */
    weislo3d slo;       /* slowness */
    weiop3d weop;       /* WEI operator */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

#ifdef _OPENMP
    ompnth=omp_init(); /* OMP parameters */
#endif

    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */
    if (!sf_getbool(  "adj",&adj   ))   adj = true; /* adjoint flag, true for migration, false for modeling */

    if(verb) fprintf(stderr,"init cube...");
    cub = wei_cube(verb,ompnth);
    if(verb) fprintf(stderr,"OK\n");

    /*------------------------------------------------------------*/
    Fslo = sf_input ("slo");
    if (SF_FLOAT !=sf_gettype(Fslo)) sf_error("Need float slowness");
    weislo_inp(cub,Fslo); /* input slowness */

    /*------------------------------------------------------------*/
    if(adj){
        Fdat = sf_input ( "in");
        if (SF_COMPLEX !=sf_gettype(Fdat)) sf_error("Need complex sdat");
    
        weidat_inp(cub,Fdat); /* input data dimensions */

        Fcic = sf_output("out"); sf_settype(Fcic,SF_FLOAT);
        weicic_out(cub,Fcic); /* prepare output image */
    }
    else{
        Fcic = sf_input("in");
        Fsou = sf_input("sou");
        if (SF_COMPLEX !=sf_gettype(Fcic)) sf_error("Need complex sdat");
        weidat_inp(cub,Fsou); /* input data dimensions */

        Fdat = sf_output("out"); sf_settype(Fdat,SF_COMPLEX);
        weidat_out(cub,Fdat); /* prepare output data */
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"init slo...");
    slo = weizoslo_init(cub,Fslo);
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
    weop = weicic_init(cub);
    if(adj)
      weizocic(weop,cub,ssr,tap,slo,Fdat,Fcic);
    else
      weizomod(weop,cub,ssr,tap,slo,Fdat,Fcic);

    weicic_close(weop);
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
