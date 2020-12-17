/* 3-D S/R WEMVA with extended split-step */
/*
  Copyright (C) 2009 Colorado School of Mines
  
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
#include "omputil.h"
#endif

#include "wex.h"
#include "wexmva.h"
#include "wexeic.h"
#include "wextap.h"
#include "wexlsr.h"
#include "wexssr.h"
#include "wexslo.h"
#include "wexutl.h"
/*^*/

int main (int argc, char *argv[])
{
    int  adj;             /* forward or adjoint */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */

    int nhx, nhy, nhz, nht, nc;
    int nhx2,nhy2,nhz2,nht2;
    float dht, oht;

    sf_axis alx,aly,az;
    sf_axis aw,ae;

    /* I/O files */
    sf_file Bwr=NULL;   /*  background wavefield file Bwr */
    sf_file Bs=NULL;    /*  background slowness file Bs   */
    sf_file Ps=NULL;    /*  slowness perturbation file Ps */
    sf_file Pi=NULL;    /*  image perturbation file Pi    */
    sf_file Pwr=NULL;   /*  perturbed wavefield file Pwr */

    int ompnth=1;

    wexcub3d cub; /* wavefield hypercube */
    wexcip3d cip; /* CIP gathers */
    wextap3d tap; /* tapering */
    wexssr3d ssr; /* SSR operator */
    wexlsr3d lsr; /* LSR operator */
    wexslo3d slo; /* slowness */

    wexmvaop3d mva;
    float dsmax;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
#endif

    if (!sf_getbool(  "verb",&verb ))  verb =  true; /* verbosity flag */
    if (!sf_getint(   "adj",&adj   ))  sf_error("Specify adjoint!"); /* y=ADJ Back-projection; n=FWD Forward Scattering */
    if (!sf_getfloat(  "eps",&eps  ))   eps =  0.01; /* stability parameter */
    if (!sf_getint(  "nrmax",&nrmax)) nrmax =     1; /* max number of refs */
    if (!sf_getfloat("dtmax",&dtmax)) dtmax = 0.004; /* max time error */
    if (!sf_getint(    "pmx",&pmx  ))   pmx =     0; /* padding on x */
    if (!sf_getint(    "pmy",&pmy  ))   pmy =     0; /* padding on y */
    if (!sf_getint(    "tmx",&tmx  ))   tmx =     0; /* taper on x   */
    if (!sf_getint(    "tmy",&tmy  ))   tmy =     0; /* taper on y   */

    ae  = sf_maxa(1,0,1);
    nhx=nhy=nhz=nht=nc=nhx2=nhy2=nhz2=nht2=0;
    oht = dht = 0.0;
    /*------------------------------------------------------------*/
    /* slowness */

    Bs = sf_input("slo");
    alx = sf_iaxa(Bs,1); sf_setlabel(alx,"lx");
    aly = sf_iaxa(Bs,2); sf_setlabel(aly,"ly");
    az =  sf_iaxa(Bs,3); sf_setlabel(az, "z");

    /*------------------------------------------------------------*/
    /* input file */
    if(adj)
        Pi = sf_input("in");
    else
        Ps = sf_input("in");

    /*------------------------------------------------------------*/
    /* wavefield */
    Bwr = sf_input("wfl");
    aw  = sf_iaxa(Bwr,4); sf_setlabel(aw ,"w" );
    Pwr = sf_tmpfile(NULL); sf_settype(Pwr,SF_COMPLEX);

    /*------------------------------------------------------------*/
    cub = wex_cube(verb,
                   alx,aly,az,
                   alx,aly,
                   aw,
                   ae,
                   eps,
                   ompnth);
    dsmax = dtmax/cub->az.d;

    /*------------------------------------------------------------*/
    /* init structures */
    tap = wextap_init(cub->alx.n,
                      cub->aly.n,
                      1,
                      SF_MIN(tmx,cub->alx.n-1), /* tmx */
                      SF_MIN(tmy,cub->aly.n-1), /* tmy */
                      0,
                      true,true,false);
    slo = wexslo_init(cub,Bs,nrmax,dsmax);
    ssr = wexssr_init(cub,slo,pmx,pmy,tmx,tmy,dsmax);
    lsr = wexlsr_init(cub,pmx,pmy,dsmax);

    /*------------------------------------------------------------*/
    /* WEMVA */
    if(adj) {
        sf_warning("adjoint operator...");

        cip = wexcip_init(cub,nhx,nhy,nhz,nht,nhx2,nhy2,nhz2,nht2,nc,dht,oht,Pi,0);
        mva = wexmva_init(cub,cip);

        Ps = sf_output("out"); sf_settype(Ps,SF_COMPLEX);
        sf_oaxa(Ps,alx,1);
        sf_oaxa(Ps,aly,2);
        sf_oaxa(Ps,az, 3);

        /* Adjoint I.C. operator, dI -> dW */
        wexzocip_adj(cub,cip,Pwr,Pi); /* Wr dR */
        sf_filefresh(Pwr);

        /* Adjoint WEMVA operator, dW -> dS */
        wexzomva(mva,adj,cub,ssr,lsr,tap,slo,Bwr,Pwr,Ps);

    } else {
        /* set up the I/O of output CIP gathers */
        Pi = sf_output("out"); sf_settype(Pi,SF_COMPLEX);

        sf_oaxa(Pi,alx,1);
        sf_oaxa(Pi,aly,2);
        sf_oaxa(Pi,az, 3);

        cip = wexcip_init(cub,nhx,nhy,nhz,nht,nhx2,nhy2,nhz2,nht2,nc,dht,oht,Pi,0);
        mva = wexmva_init(cub,cip);
 
        /* WEMVA operator, dS -> dW */
        wexzomva(mva,adj,cub,ssr,lsr,tap,slo,Bwr,Pwr,Ps);
        sf_filefresh(Pwr);

        /* I.C. operator, dW -> dI */
        wexzocip_for(cub,cip,Pwr,Pi); /* CONJ( Ws) dWr */

    }

    /*------------------------------------------------------------*/
    /* close structures */
    wexslo_close(slo);
    wexssr_close(cub,ssr);
    wextap2D_close(tap);
    wexmva_close(mva);
    wexcip_close(cip,0);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* close files */
    if (Ps!=NULL) sf_fileclose(Ps);
    if (Pi!=NULL) sf_fileclose(Pi);
    if (Bwr!=NULL) sf_fileclose(Bwr);
    if (Pwr!=NULL) sf_tmpfileclose(Pwr);
    /*------------------------------------------------------------*/

    exit (0);
}

