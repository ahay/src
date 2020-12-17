/* 3-D zero-offset modeling/migration with extended SSF */

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
#include "wexeic.h"
#include "wextap.h"
#include "wexssr.h"
#include "wexslo.h"
#include "wexutl.h"
/*^*/

int main (int argc, char *argv[])
{
    bool verb;            /* verbosity */
    bool causal;          /* causal/anti-causal flag */
    float eps;            /* dip filter constant */
    int  adj;             /* forward or adjoint */
    int save;             /* save wavefield or not */
    int   nrmax;          /* number of reference velocities */
    float dtmax;          /* time error */
    int   pmx,pmy;        /* padding in the k domain */
    int   tmx,tmy;        /* boundary taper size */
    int wexsign;

    int nw, iw;

    int nhx, nhy, nhz, nht, nc;
    int nhx2,nhy2,nhz2,nht2;
    float dht, oht, ow, dw;

    sf_axis alx,aly,az;
    sf_axis aw,ae;

    /* I/O files */
    sf_file Fs=NULL;    /*  slowness file S(nlx,nly,nz   ) */            
    sf_file Fd=NULL;    /*  if adj=0, then modeling, Fd is 
			    output surface data file D(nmx,nmy,nw); 
			    otherwise Fd is input surface data 
                            file D(nmx,nmy,nw)             */
    sf_file Fm=NULL;    /*  if adj=0, Fm is input reflectivity 
                            r(nmx,nmy,nmz); otherwise Fm is output 
			    image r(nmx,nmy,nmz) or CIP 
			    gathers                        */
    sf_file Fwr=NULL;   /*  receiver wavefield             */
    sf_file Fis=NULL;   /*  injected source wavefield      */
    sf_file Ftw=NULL;   /*  temp wavefield                 */

    /* temp arrays */
    sf_complex ***wtmp;

    int ompnth=1;

    wexcub3d cub; /* wavefield hypercube */
    wexcip3d cip; /* CIP gathers */
    wextap3d tap; /* tapering */
    wexssr3d ssr; /* SSR operator */
    wexslo3d slo; /* slowness */

    wexop3d weop;

    float dsmax;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
#endif

    if (!sf_getint(    "adj",&adj   ))  sf_error("Specify the adjoint!"); /* y=ADJ migration; n=FWD modeling */
    if (!sf_getint(   "save",&save  ))  save = 0;     /* save wfld flag */
    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */
    if (!sf_getfloat(  "eps",&eps   ))   eps =  0.01; /* stability parameter */
    if (!sf_getint(  "nrmax",&nrmax )) nrmax =     1; /* maximum references */
    if (!sf_getfloat("dtmax",&dtmax )) dtmax = 0.004; /* max time error */

    if (!sf_getint(    "pmx",&pmx   ))   pmx =     0; /* padding on x */
    if (!sf_getint(    "pmy",&pmy   ))   pmy =     0; /* padding on y*/

    if (!sf_getint(    "tmx",&tmx   ))   tmx =     0; /* taper on x*/
    if (!sf_getint(    "tmy",&tmy   ))   tmy =     0; /* taper on y */

    if(adj) causal= false; else causal= true;

    if(!adj){
        Fm = sf_input("in");
        if (!sf_getint  ("nw",&nw)) sf_error ("Need nw=");
        if (!sf_getfloat("dw",&dw)) sf_error ("Need dw=");
        if (!sf_getfloat("ow",&ow)) ow=0.;
        aw = sf_maxa(nw,ow,dw);
        sf_setunit (aw,"1/s");
        sf_setlabel(aw , "w");
    }
    else{
        Fd = sf_input(  "in");
        aw  = sf_iaxa(Fd,3); 
    }

    ae  = sf_maxa(1,0,1);
    nhx=nhy=nhz=nht=nc=nhx2=nhy2=nhz2=nht2=0;
    oht = dht = 0.0;

    /*------------------------------------------------------------*/
    /* slowness */

    Fs = sf_input("slo");
    alx = sf_iaxa(Fs,1); sf_setlabel(alx,"x");
    aly = sf_iaxa(Fs,2); sf_setlabel(aly,"y");
    az =  sf_iaxa(Fs,3); sf_setlabel(az, "z"); sf_setunit(az, "km");

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
    slo = wexslo_init(cub,Fs,nrmax,dsmax);
    ssr = wexssr_init(cub,slo,pmx,pmy,tmx,tmy,dsmax);
    weop = wex_init(cub);

    wexsign=causal?+1:-1;

    /*------------------------------------------------------------*/
    Ftw = sf_tmpfile(NULL); sf_settype(Ftw,SF_COMPLEX);
    Fis = sf_tmpfile(NULL); sf_settype(Fis,SF_COMPLEX);

    /*------------------------------------------------------------*/
    if(!adj){ /* For modeling, output surface data */
        sf_warning("modeling ....");

        cip = wexcip_init(cub,nhx,nhy,nhz,nht,nhx2,nhy2,nhz2,nht2,nc,dht,oht,Fd,0);

	Fd = sf_output("out"); sf_settype(Fd,SF_COMPLEX);
        sf_oaxa(Fd,alx,1);
        sf_oaxa(Fd,aly,2);
        sf_oaxa(Fd,aw, 3);

	wexzocip_adj(cub,cip,Fis,Fm);
        sf_filefresh(Fis);

        // wavefield extrapolation 
        wex(weop,cub,ssr,tap,slo,wexsign,Fis,Ftw,adj);

        // allocation for temp arrays 
        wtmp = sf_complexalloc3(cub->amx.n,cub->aly.n,cub->aw.n);

        for(iw=0; iw<cub->aw.n; iw++){
            sf_seek(Ftw,sizeof(sf_complex)*cub->amx.n*cub->aly.n*cub->az.n*iw,SEEK_SET);
            sf_complexread(wtmp[iw][0],cub->amx.n*cub->aly.n,Ftw);
        }
        sf_complexwrite(wtmp[0][0],cub->amx.n*cub->aly.n*cub->aw.n,Fd);

        /* deallocate temp arrays */
        free(**wtmp );free(*wtmp );free( wtmp );
    }
    else{ /* For migration, output image or CIPs */
        sf_warning("migration ....");

        if (SF_COMPLEX !=sf_gettype(Fd)) sf_error("Need complex data");
        Fm = sf_output("out"); sf_settype(Fm,SF_COMPLEX);

	sf_oaxa(Fm,alx,1);
	sf_oaxa(Fm,aly,2);
	sf_oaxa(Fm,az, 3);

        /* zero-offset wavefield extrapolation */
        wex(weop,cub,ssr,tap,slo,wexsign,Fd,Ftw,adj);

        /* initialize image gathers */
        cip = wexcip_init(cub,nhx,nhy,nhz,nht,nhx2,nhy2,nhz2,nht2,nc,dht,oht,Fd,0);

        /* construct image */
        wexzocip_for(cub,cip,Ftw,Fm);

        if(save){
            Fwr = sf_output("wfl"); sf_settype(Fwr,SF_COMPLEX);
            sf_oaxa(Fwr,alx,1);
            sf_oaxa(Fwr,aly,2);
            sf_oaxa(Fwr,az, 3);
            sf_oaxa(Fwr,aw, 4);

            sf_filefresh(Ftw);
            sf_filecopy(Fwr,Ftw,SF_COMPLEX);
        }
    }

    /*------------------------------------------------------------*/
    /* close structures */
    wex_close(weop);
    wexslo_close(slo);
    wexssr_close(cub,ssr);
    wextap2D_close(tap);
    wexcip_close(cip,0);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* close files */
    if (Fs!=NULL) sf_fileclose(Fs);
    if (Fd!=NULL) sf_fileclose(Fd);
    if (Fm!=NULL) sf_fileclose(Fm);
    if (Fwr!=NULL) sf_fileclose(Fwr);
    if (Fis!=NULL) sf_tmpfileclose(Fis);
    if (Ftw!=NULL) sf_tmpfileclose(Ftw);
    /*------------------------------------------------------------*/

    exit (0);
}
         



