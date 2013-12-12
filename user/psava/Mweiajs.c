/* Adjoint source construction for image-domain waveform tomography */

/*
  Copyright (C) 2011 Colorado School of Mines

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
#include "weiutl.h"
/*^*/

#include "weiajs.h"

int main (int argc, char *argv[])
{
    bool  verb;         /* verbosity */
    int ompnth=1;       /* number of threads */
    bool conj;          /* complex conjugate flag */
    char *irun;         /* I run... */

    sf_file Faso=NULL;  /* output     wavefield wfl(x,y,z,w) */
    sf_file Fbwf=NULL;  /*  input     wavefield wfl(x,y,z,w) */
    sf_file Feic=NULL;  /* EIC image file   R(hx,hy,hz,tt,c) */
    sf_file Fcoo=NULL;  /*     coord file   C(nc) */

    weicub3d cub;       /* hypercube */
    weiop3d weop;       /* WEI operator */
    weico3d eico;       /* eic coordinates */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

#ifdef _OPENMP
    ompnth=omp_init(); /* OMP parameters */
#endif

    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */
    if (!sf_getbool(  "conj",&conj   )) sf_error("Specify whether complex conjugate!"); /* flag */
    if (NULL == (irun = sf_getstring("irun"))) irun = "e";

    if(verb) fprintf(stderr,"init cube...");
    cub = wei_cube(verb,ompnth);
    if(verb) fprintf(stderr,"OK\n");

    /*------------------------------------------------------------*/
    Fbwf = sf_input("bwf");
    if (SF_COMPLEX !=sf_gettype(Fbwf)) sf_error("Need complex wavefield");
    weiwfl_inp(cub,Fbwf); /* input data dimensions */

    /*------------------------------------------------------------*/
    Feic = sf_input( "in");
  
    /*------------------------------------------------------------*/
    Fcoo = sf_input ("coo"); /* input coordinates for EIC */
    if (SF_FLOAT !=sf_gettype(Fcoo)) sf_error("Need float coordinates");

    /* load EIC dimensions from file */
    mvahic_inp(cub,Feic,Fcoo); 

    eico = weicoo_init(cub,Fcoo);

    /*------------------------------------------------------------*/
    Faso = sf_output("out"); sf_settype(Faso,SF_COMPLEX);
    weiwfl_out(cub,Faso);

    /*------------------------------------------------------------*/
    switch(irun[0]) {

	case 'h': /* HIC: hx-hy-1  */
	    if(verb) sf_warning("HIC AJS");

	    weop = weihicajs_init(cub);
 
	    weihicajs(weop,cub,eico,Feic,Fbwf,Faso,conj);

	    weihicajs_close(weop);
	    
	    break;

	case 'e': /* EIC  hx-hy-hz */
	default:
	    if(verb) sf_warning("EIC AJS");

	    adjsou(cub,eico,Feic,Fbwf,Faso,conj);

	    break;
    }
	
    weicoo_close(eico);
    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* close files */
    if (Fbwf!=NULL) sf_fileclose(Fbwf);
    if (Faso!=NULL) sf_fileclose(Faso);
    if (Fcoo!=NULL) sf_fileclose(Fcoo);
    if (Feic!=NULL) sf_fileclose(Feic);
    /*------------------------------------------------------------*/

    exit (0);
}

