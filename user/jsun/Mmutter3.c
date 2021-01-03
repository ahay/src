/* 3D muting of scalor or vector data */
/*
  Copyright (C) 2016 The University of Texas at Austin
  
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

  The main framework is following Paul Sava's sfewefd3d program.
  This program uses pseudo-spectral method to calculate spatial derivatives.
*/

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ksutil.h"

int main(int argc, char* argv[])
{
    /*------------------------------------------------------------*/
    /* Execution control, I/O files and geometry                  */
    /*------------------------------------------------------------*/
    bool verb; /* execution flags */

    /* I/O files */
    sf_file Finp=NULL; /* input     */
    sf_file Fout=NULL; /* output    */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */

    /* cube axes */
    sf_axis at,ax,ay; /* time, x, y, z */ 
    sf_axis asx,asy,arx,ary,ac;    /* sou, rec-x, rec-y, component */ 

    /* dimension, index and interval */
    int nt,ns,nr,nc;

    /* data and mute structure */
    dat3d dat=NULL;
    mut3d mut=NULL;

    /* I/O arrays for sou & rec */
    pt3d    *ss=NULL;           /* sources   */
    pt3d    *rr=NULL;           /* receivers */
    float ***dd=NULL;           /* data      */

    /* mutting related */
    float t0,velw,eps;

    /*------------------------------------------------------------*/
    /* init RSF                                                   */
    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters                                             */
    /*------------------------------------------------------------*/
#ifdef _OPENMP
    omp_init();
#endif

    /*------------------------------------------------------------*/
    /* read execution flags and pars                              */
    /*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getfloat("t0",&t0)) t0=0; /* source delay time */
    if(! sf_getfloat("velw",&velw)) velw=1.5; /* water velocity */
    if(! sf_getfloat("eps",&eps)) eps=1000.; /* decay parameter */

    /*------------------------------------------------------------*/
    /* I/O files                                                  */
    /*------------------------------------------------------------*/
    Finp = sf_input ("in" ); /* input  */
    Fsou = sf_input ("sou"); /* sources */
    Frec = sf_input ("rec"); /* receivers */
    Fout = sf_output("out"); /* output */

    /*------------------------------------------------------------*/
    /* axes                                                       */
    /*------------------------------------------------------------*/
    ax = sf_iaxa(Finp,1); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space x */
    ay = sf_iaxa(Finp,2); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space y */
    ac = sf_iaxa(Finp,3); sf_setlabel(ac,"c"); if(verb) sf_raxa(ac); /* component */
    at = sf_iaxa(Finp,4); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */

    asx = sf_iaxa(Fsou,2); sf_setlabel(asx,"sx"); if(verb) sf_raxa(asx); /* sources x */
    asy = sf_iaxa(Fsou,3); sf_setlabel(asy,"sy"); if(verb) sf_raxa(asy); /* sources y */
    arx = sf_iaxa(Frec,2); sf_setlabel(arx,"rx"); if(verb) sf_raxa(arx); /* receivers x */
    ary = sf_iaxa(Frec,3); sf_setlabel(ary,"ry"); if(verb) sf_raxa(ary); /* receivers y */

    nt = sf_n(at);

    nc = sf_n(ac);

    ns = sf_n(asx)*sf_n(asy);
    nr = sf_n(arx)*sf_n(ary);

    /*------------------------------------------------------------*/
    /* set up data domain parameter                               */
    /*------------------------------------------------------------*/
    dat=dat3d_init(ns,nr,nc,at,ax,ay);

    /*------------------------------------------------------------*/
    /* setup output data header                                   */
    /*------------------------------------------------------------*/
    sf_oaxa(Fout,arx,1);
    sf_oaxa(Fout,ary,2);
    sf_oaxa(Fout,ac,3);
    sf_oaxa(Fout,at,4);

    /*------------------------------------------------------------*/
    /* source and data array                                      */
    /*------------------------------------------------------------*/
    dd=sf_floatalloc3(nr,nc,nt);
    sf_floatread(dd[0][0],nr*nc*nt,Finp);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates                          */
    /*------------------------------------------------------------*/
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */

    /* calculate 3d mutting coef for sou & rec */
    mut = mut3d_make(t0,velw,eps,ss,rr,dat);

    /*------------------------------------------------------------*/ 
    /* mute and output                                            */ 
    /*------------------------------------------------------------*/
    mut3d_apply(dd,dat,mut);
    sf_floatwrite(dd[0][0],nr*nc*nt,Fout);
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    
    free(ss);
    free(rr);
    free(**dd); free(*dd); free(dd);

    /*------------------------------------------------------------*/

    exit (0);
}

