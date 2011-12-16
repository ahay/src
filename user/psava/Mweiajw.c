/* 3-D wave-equation wavefield continuation with adjoint-source */

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
#include "weissr.h"
#include "weiutl.h"
#include "weitap.h"
#include "weislo.h"

int main (int argc, char *argv[])
{
    bool  verb;         /* verbosity */
    int ompnth=1;       /* number of threads */
    bool  down;         /* adjoint flag */
    bool causal;

    sf_file Fslo=NULL;  /*  slowness file  S(nlx,nly,nz   ) */
    sf_file Fsou=NULL;  /*      data file Ds(nmx,nmy,   nw) */
    sf_file Fwfl=NULL;  /*      data file Dr(nmx,nmy,   nw) */

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
    if (!sf_getbool(  "down",&down  ))  down = true;  /* up/down   flag */
    if (!sf_getbool("causal",&causal))  sf_error("Specify causal!"); /* causality flag */

    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"init cube...");
    cub = wei_cube(verb,ompnth);
    if(verb) fprintf(stderr,"OK\n");

    /*------------------------------------------------------------*/
    Fslo = sf_input ("slo");
    if (SF_FLOAT !=sf_gettype(Fslo)) sf_error("Need float slowness");
    weislo_inp(cub,Fslo); /* input slowness */

    /*------------------------------------------------------------*/
    Fsou = sf_input ( "in");
    if (SF_COMPLEX !=sf_gettype(Fsou)) sf_error("Need complex sdat");
    weiwfl_inp(cub,Fsou); /* input adjoint source dimensions */

    /*------------------------------------------------------------*/
    Fwfl = sf_output("out"); sf_settype(Fwfl,SF_COMPLEX);

    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"init slo...");
    slo = weislo_init(cub,Fslo);
    if(verb) fprintf(stderr,"OK\n");
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
    if(verb) fprintf(stderr,"init weop...");
    weop = weiwfl_init(cub);
    if(verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    adjwfl(weop,cub,ssr,tap,slo,Fsou,Fwfl,down,causal);
    weiwfl_close(weop);

    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"close structures...");
    weislo_close(slo);
    weissr_close(cub,ssr);
    weitap_close(tap);
    if(verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* close files */
    if (Fsou!=NULL) sf_fileclose(Fsou);
    if (Fwfl!=NULL) sf_fileclose(Fwfl);
    if (Fslo!=NULL) sf_fileclose(Fslo);
    /*------------------------------------------------------------*/

    exit (0);
}

