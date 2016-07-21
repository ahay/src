/* 3D elastic recursive integral time extrapolation using KISS-FFT
   sou wavelet  (nx,ny,nc,nt)
   rec data     (nx,ny,nc,nt)
   sou geometry (nc,nx,ny)
   rec geometry (nc,nx,ny)
*/
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
    sf_file Finp=NULL; /* wavefield */
    sf_file Fout=NULL; /* wavefield */
    sf_file Frnk=NULL; /* app. rank */
    sf_file Flft=NULL; /* left mat  */
    sf_file Frht=NULL; /* right mat */

    /* cube axes */
    sf_axis ax,ay,az,ac; /* time, x, y, z */ 

    /* dimension, index and interval */
    int     nz,nx,ny,nb,nc;
    float   dz,dx,dy;
    int     nxyz, nk;

    /* FDM and KSP structure */ //!!!JS
    fdm3d    fdm=NULL;
    dft3d    dft=NULL;
    clr3d    clr=NULL;

    /*------------------------------------------------------------*/
    /* displacement: uo = U @ t; up = U @ t+1                     */
    /*------------------------------------------------------------*/
    sf_complex ***uox, ***uoy, ***uoz, **uo;

    /*------------------------------------------------------------*/
    /* lowrank decomposition arrays                               */
    /*------------------------------------------------------------*/
    int ntmp, *n2s;
    sf_complex **lt, **rt;

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
    /* read execution flags                                       */
    /*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */

    /*------------------------------------------------------------*/
    /* I/O files                                                  */
    /*------------------------------------------------------------*/
    Finp = sf_input ("in" ); /* input     */
    Frnk = sf_input ("rnk"); /* app. rank */
    Flft = sf_input ("lft"); /* left mat  */
    Frht = sf_input ("rht"); /* right mat */
    Fout = sf_output("out"); /* output    */

    /*------------------------------------------------------------*/
    /* axes                                                       */
    /*------------------------------------------------------------*/
    az = sf_iaxa(Finp,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Finp,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space x */
    ay = sf_iaxa(Finp,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space y */

    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);
    nb = 0;

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC                     */
    /*------------------------------------------------------------*/
    fdm=fdutil3d_init(verb,false,az,ax,ay,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);

    /*------------------------------------------------------------*/
    /* 3D vector components                                       */
    /*------------------------------------------------------------*/
    nc=3;
    ac=sf_maxa(nc,0,1); /* output 3 cartesian components */

    /*------------------------------------------------------------*/
    /* setup output header                                        */
    /*------------------------------------------------------------*/
    sf_settype(Fout,SF_COMPLEX);
    sf_oaxa(Fout,az,1);
    sf_oaxa(Fout,ax,2);
    sf_oaxa(Fout,ay,3);
    sf_oaxa(Fout,ac,4);

    /*------------------------------------------------------------*/
    /* allocate and initialize wavefield arrays                   */
    /*------------------------------------------------------------*/
    /* z-component */
    uoz=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    /* x-component */
    uox=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    /* y-component */
    uoy=sf_complexalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    /* wavefield vector */
    uo = (sf_complex**) sf_alloc(3,sizeof(sf_complex*));
    uo[0] = uox[0][0]; uo[1] = uoy[0][0]; uo[2] = uoz[0][0];

    /* initialize fft and lrk */
    dft = dft3d_init(1,false,false,fdm);
    nxyz= fdm->nypad*fdm->nxpad*fdm->nzpad;
    nk  = dft->nky  *dft->nkx  *dft->nkz;

    /*------------------------------------------------------------*/ 
    /* allocation I/O arrays                                      */
    /*------------------------------------------------------------*/ 
    n2s = sf_intalloc(6);
    sf_intread(n2s,6,Frnk);
    clr = clr3d_make(n2s,fdm);
    clr3d_init(fdm,dft,clr);

    /* check the dimension */
    if (!sf_histint(Flft,"n1",&ntmp) || ntmp != nxyz)        sf_error("Need n1=%d in left",nxyz);
    if (!sf_histint(Flft,"n2",&ntmp) || ntmp != clr->n2_sum) sf_error("Need n2=%d in left",clr->n2_sum);
    if (!sf_histint(Frht,"n1",&ntmp) || ntmp != nk)          sf_error("Need n1=%d in right",nk);
    if (!sf_histint(Frht,"n2",&ntmp) || ntmp != clr->n2_sum) sf_error("Need n2=%d in right",clr->n2_sum);
  
    lt  = sf_complexalloc2(nxyz,clr->n2_sum); 
    rt  = sf_complexalloc2(nk  ,clr->n2_sum); 
    sf_complexread(lt[0],nxyz*clr->n2_sum,Flft);
    sf_complexread(rt[0],nk  *clr->n2_sum,Frht);

    sf_complexread(uoz[0][0],nxyz,Finp);
    sf_complexread(uox[0][0],nxyz,Finp);
    sf_complexread(uoy[0][0],nxyz,Finp);

    /*------------------------------------------------------------*/
    /* apply lowrank matrix to wavefield vector                   */
    /*------------------------------------------------------------*/
    clr3d_apply(uo, uo, lt, rt, fdm, dft, clr);

    /*------------------------------------------------------------*/
    /* cut wavefield and save */
    /*------------------------------------------------------------*/
    sf_complexwrite(uoz[0][0],nxyz,Fout);
    sf_complexwrite(uox[0][0],nxyz,Fout);
    sf_complexwrite(uoy[0][0],nxyz,Fout);

    /*------------------------------------------------------------*/
    /* deallocate arrays */

    free(dft);
    dft3d_finalize();
    free(clr);
    clr3d_finalize();

    free(n2s);
    free(*lt); free(lt); free(*rt); free(rt);

    free(**uoz); free(*uoz); free(uoz);
    free(uo);

    /*------------------------------------------------------------*/

    exit (0);
}

