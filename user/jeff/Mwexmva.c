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
    int  eic;             /* EIC or CIC */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant */
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */

    int nhx, nhy, nhz, nht, nc;
    int nhx2,nhy2,nhz2,nht2;
    float dht, oht;

    sf_axis amx,amy,az;
    sf_axis alx,aly;
    sf_axis aw,ae,ac,aa;
    sf_axis ahx,ahy,ahz,aht;         

    /* I/O files */
    sf_file Bws=NULL;   /*  background wavefield file Bws */
    sf_file Bwr=NULL;   /*  background wavefield file Bwr */
    sf_file Bs=NULL;    /*  background slowness file Bs   */
    sf_file Ps=NULL;    /*  slowness perturbation file Ps */
    sf_file Pi=NULL;    /*  image perturbation file Pi    */
    sf_file Fc=NULL;    /*  CIP coordinates               */
    sf_file Pws=NULL;   /*  perturbed wavefield file Pws */
    sf_file Pwr=NULL;   /*  perturbed wavefield file Pwr */
    sf_file Pti=NULL;  

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
    if (!sf_getint(  "feic",&eic   ))  sf_error("Specify EIC!");     /* extended I.C. flag */
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
    Bws = sf_input("swfl");
    Bwr = sf_input("rwfl");

    amx = sf_iaxa(Bws,1); sf_setlabel(amx,"mx");
    amy = sf_iaxa(Bws,2); sf_setlabel(amy,"my");
    aw  = sf_iaxa(Bws,4); sf_setlabel(aw ,"w" );

    Pws = sf_tmpfile(NULL); sf_settype(Pws,SF_COMPLEX);
    Pwr = sf_tmpfile(NULL); sf_settype(Pwr,SF_COMPLEX);

    /*------------------------------------------------------------*/
    cub = wex_cube(verb,
                   amx,amy,az,
                   alx,aly,
                   aw,
                   ae,
                   eps,
                   ompnth);
    dsmax = dtmax/cub->az.d;

    /*------------------------------------------------------------*/
    /* init structures */
    tap = wextap_init(cub->amx.n,
                      cub->amy.n,
                      1,
                      SF_MIN(tmx,cub->amx.n-1), /* tmx */
                      SF_MIN(tmy,cub->amy.n-1), /* tmy */
                      0,
                      true,true,false);
    slo = wexslo_init(cub,Bs,nrmax,dsmax);
    ssr = wexssr_init(cub,slo,pmx,pmy,tmx,tmy,dsmax);
    lsr = wexlsr_init(cub,pmx,pmy,dsmax);

    /*------------------------------------------------------------*/
    Pti = sf_tmpfile(NULL); sf_settype(Pti,SF_COMPLEX);

    /*------------------------------------------------------------*/
    /* WEMVA */
    if(adj) {
        sf_warning("adjoint operator...");

        if(eic){
            ahx = sf_iaxa(Pi,1); sf_setlabel(ahx,"hx");
            ahy = sf_iaxa(Pi,2); sf_setlabel(ahy,"hy");
            ahz = sf_iaxa(Pi,3); sf_setlabel(ahz,"hz");
            aht = sf_iaxa(Pi,4); sf_setlabel(aht,"ht");

            dht = sf_d(aht);  oht = sf_o(aht);

            nhx2 = sf_n(ahx); nhx = (nhx2-1)/2;
            nhy2 = sf_n(ahy); nhy = (nhy2-1)/2;
            nhz2 = sf_n(ahz); nhz = (nhz2-1)/2;
            nht2 = sf_n(aht); nht = (nht2-1)/2;

            /* CIP coordinates */
            Fc = sf_input ("cc" );
            ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
            nc = sf_n(ac);
        }

        cip = wexcip_init(cub,nhx,nhy,nhz,nht,nhx2,nhy2,nhz2,nht2,nc,dht,oht,Fc,eic);
        mva = wexmva_init(cub,cip);

        Ps = sf_output("out"); sf_settype(Ps,SF_COMPLEX);
        sf_oaxa(Ps,amx,1);
        sf_oaxa(Ps,amy,2);
        sf_oaxa(Ps,az, 3);
        if(eic){
        sf_oaxa(Ps,ae, 4);
        sf_oaxa(Ps,ae, 5);}

        /* Adjoint I.C. operator, dI -> dW */
        wexcip_adj(cub,cip,Bwr,Pws,Pi,eic,1,1); /* Ws dR */
        wexcip_adj(cub,cip,Bws,Pwr,Pi,eic,0,0); /* Wr dR */

        sf_filefresh(Pws);
        sf_filefresh(Pwr);

        /* Adjoint WEMVA operator, dW -> dS */
        wexmva(mva,adj,cub,ssr,lsr,tap,slo,Bws,Bwr,Pws,Pwr,Ps);

    } else {
        /* set up the I/O of output CIP gathers */
        Pi = sf_output("out"); sf_settype(Pi,SF_COMPLEX);

        if(eic){
            /* CIP coordinates */
            Fc = sf_input ("cc" );
            ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
            nc = sf_n(ac);

            if(! sf_getint("nhx",&nhx)) nhx=0; /* number of lags on the x axis */
            if(! sf_getint("nhy",&nhy)) nhy=0; /* number of lags on the y axis */
            if(! sf_getint("nhz",&nhz)) nhz=0; /* number of lags on the z axis */
            if(! sf_getint("nht",&nht)) nht=0; /* number of lags on the t axis */
            if(! sf_getfloat("dht",&dht)) sf_error("need dht");
            oht = -nht*dht;

            nhx2=2*nhx+1; nhy2=2*nhy+1;
            nhz2=2*nhz+1; nht2=2*nht+1;

            aa=sf_maxa(nhx2,-nhx*cub->amx.d,cub->amx.d);
            sf_setlabel(aa,"hx"); sf_setunit(aa,"");
            if(verb) sf_raxa(aa);
            sf_oaxa(Pi,aa,1);

            aa=sf_maxa(nhy2,-nhy*cub->amy.d,cub->amy.d);
            sf_setlabel(aa,"hy"); sf_setunit(aa,"");
            if(verb) sf_raxa(aa);
            sf_oaxa(Pi,aa,2);

            aa=sf_maxa(nhz2,-nhz*cub->az.d,cub->az.d);
            sf_setlabel(aa,"hz"); sf_setunit(aa,"");
            if(verb) sf_raxa(aa);
            sf_oaxa(Pi,aa,3);

            aa=sf_maxa(nht2,-nht*dht,dht);
            sf_setlabel(aa,"ht"); sf_setunit(aa,"s");
            if(verb) sf_raxa(aa);
            sf_oaxa(Pi,aa,4);
 
            sf_oaxa(Pi,ac,5);

        }
        else{
            sf_oaxa(Pi,amx,1);
            sf_oaxa(Pi,amy,2);
            sf_oaxa(Pi,az, 3);
        }

        cip = wexcip_init(cub,nhx,nhy,nhz,nht,nhx2,nhy2,nhz2,nht2,nc,dht,oht,Fc,eic);
        mva = wexmva_init(cub,cip);
 
        /* WEMVA operator, dS -> dW */
        wexmva(mva,adj,cub,ssr,lsr,tap,slo,Bws,Bwr,Pws,Pwr,Ps);
         
        sf_filefresh(Pws);
        sf_filefresh(Pwr);

        /* I.C. operator, dW -> dI */
        wexcip_for(cub,cip,Bws,Pwr,Pti,eic,0,0); /* CONJ( Ws) dWr */
        sf_seek(Pti,(off_t)0,SEEK_SET);
        wexcip_for(cub,cip,Pws,Bwr,Pti,eic,0,1); /* CONJ(dWs)  Wr */

        sf_filefresh(Pti);
        sf_filecopy(Pi,Pti,SF_COMPLEX);
    }

    /*------------------------------------------------------------*/
    /* close structures */
    wexslo_close(slo);
    wexssr_close(cub,ssr);
    wextap2D_close(tap);
    wexmva_close(mva);
    wexcip_close(cip,eic);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* close files */
    if (Ps!=NULL) sf_fileclose(Ps);
    if (Fc!=NULL) sf_fileclose(Fc);
    if (Pi!=NULL) sf_fileclose(Pi);
    if (Bws!=NULL) sf_fileclose(Bws);
    if (Bwr!=NULL) sf_fileclose(Bwr);
    if (Pws!=NULL) sf_tmpfileclose(Pws);
    if (Pwr!=NULL) sf_tmpfileclose(Pwr);
    if (Pti!=NULL) sf_tmpfileclose(Pti);
    /*------------------------------------------------------------*/

    exit (0);
}

