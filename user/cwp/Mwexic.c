/* Imaging Condition for WEXMIG */

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

#include "wexeic.h"
#include "wexutl.h"
/*^*/

int main (int argc, char *argv[])
{

    bool verb;            /* verbosity */
    bool  all;            /* output both ic or not */
    bool  eic;            /* output cic or eic cube */
    bool  drv;            /* output derivative */

    int nhx, nhy, nhz, nht, nc;
    int nhx2,nhy2,nhz2,nht2;
    float dht, oht;

    sf_axis amx,amy,az;
    sf_axis aw,ae,aa,ac;

    /* I/O files */
    sf_file Fws=NULL;   /*  source wavefield                */
    sf_file Fwr=NULL;   /*  receiver wavefield              */
    sf_file Fc=NULL;
    sf_file Fi=NULL;
    sf_file Fe=NULL;

    int ompnth=1;

    wexcub3d cub; /* wavefield hypercube */
    wexcip3d cip; /* CIP gathers */


    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
#endif

    if (!sf_getbool(  "fall",&all   ))   sf_error("Specify output both IC or not!");
    if (!sf_getbool(  "feic",&eic   ))   sf_error("Specify output eIC or cIC!");
    if (!sf_getbool(  "fdrv",&drv   ))   drv = false;
    if (!sf_getbool(  "verb",&verb  ))  verb = false; /* verbosity flag */

    ae  = sf_maxa(1,0,1);
    nhx=nhy=nhz=nht=nc=nhx2=nhy2=nhz2=nht2=0;
    oht = 0.0;

    /*------------------------------------------------------------*/
    /* wavefields */
    Fws = sf_input ( "in");
    Fwr = sf_input ( "rwfl");

    amx = sf_iaxa(Fwr,1); sf_setlabel(amx,"mx");
    amy = sf_iaxa(Fwr,2); sf_setlabel(amy,"my");
    az =  sf_iaxa(Fwr,3); sf_setlabel(az, "z");
    aw  = sf_iaxa(Fwr,4); sf_setlabel(aw ,"w" );

    /*------------------------------------------------------------*/
    cub = wex_cube(verb,
                   amx,amy,az,
                   amx,amy,
                   aw,
                   ae,
                   0.00005,
                   ompnth);

    /*------------------------------------------------------------*/
    /* output */
    if(all){
      Fi = sf_output("out"); sf_settype(Fi,SF_COMPLEX);
      Fe = sf_output("cip"); sf_settype(Fe,SF_COMPLEX);

      sf_oaxa(Fi,amx,1);
      sf_oaxa(Fi,amy,2);
      sf_oaxa(Fi,az, 3);
      sf_oaxa(Fi,ae, 4);

      /* CIP coordinates */
      Fc = sf_input ("cc" );
      ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
      nc = sf_n(ac);

      Fe = sf_output("cip"); sf_settype(Fe,SF_COMPLEX);

      if(! sf_getint("nhx",&nhx)) nhx=0; /* number of lags on the x axis */
      if(! sf_getint("nhy",&nhy)) nhy=0; /* number of lags on the y axis */
      if(! sf_getint("nhz",&nhz)) nhz=0; /* number of lags on the z axis */
      if(! sf_getint("nht",&nht)) nht=0; /* number of lags on the t axis */
      if(! sf_getfloat("dht",&dht)) dht=0.004;;
      oht = -nht*dht;

      nhx2=2*nhx+1;
      nhy2=2*nhy+1;
      nhz2=2*nhz+1;
      nht2=2*nht+1;

      aa=sf_maxa(nhx2,-nhx*cub->amx.d,cub->amx.d);
      sf_setlabel(aa,"hx"); sf_setunit(aa,"");
      if(verb) sf_raxa(aa);
      sf_oaxa(Fe,aa,1);

      aa=sf_maxa(nhy2,-nhy*cub->amy.d,cub->amy.d);
      sf_setlabel(aa,"hy"); sf_setunit(aa,"");
      if(verb) sf_raxa(aa);
      sf_oaxa(Fe,aa,2);

      aa=sf_maxa(nhz2,-nhz*cub->az.d,cub->az.d);
      sf_setlabel(aa,"hz"); sf_setunit(aa,"");
      if(verb) sf_raxa(aa);
      sf_oaxa(Fe,aa,3);

      aa=sf_maxa(nht2,-nht*dht,dht);
      sf_setlabel(aa,"ht"); sf_setunit(aa,"s");
      if(verb) sf_raxa(aa);
      sf_oaxa(Fe,aa,4);

      sf_oaxa(Fe,ac,5);
    }
    else{
      if(eic){
        /* CIP coordinates */
        Fc = sf_input ("cc" );
        ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
        nc = sf_n(ac);

        Fe = sf_output("out"); sf_settype(Fe,SF_COMPLEX);

        if(! sf_getint("nhx",&nhx)) nhx=0; /* number of lags on the x axis */
        if(! sf_getint("nhy",&nhy)) nhy=0; /* number of lags on the y axis */
        if(! sf_getint("nhz",&nhz)) nhz=0; /* number of lags on the z axis */
        if(! sf_getint("nht",&nht)) nht=0; /* number of lags on the t axis */
        if(! sf_getfloat("dht",&dht)) dht=0.004;;
        oht = -nht*dht;

        nhx2=2*nhx+1;
        nhy2=2*nhy+1;
        nhz2=2*nhz+1;
        nht2=2*nht+1;

        aa=sf_maxa(nhx2,-nhx*cub->amx.d,cub->amx.d);
        sf_setlabel(aa,"hx"); sf_setunit(aa,"");
        if(verb) sf_raxa(aa);
        sf_oaxa(Fe,aa,1);

        aa=sf_maxa(nhy2,-nhy*cub->amy.d,cub->amy.d);
        sf_setlabel(aa,"hy"); sf_setunit(aa,"");
        if(verb) sf_raxa(aa);
        sf_oaxa(Fe,aa,2);

        aa=sf_maxa(nhz2,-nhz*cub->az.d,cub->az.d);
        sf_setlabel(aa,"hz"); sf_setunit(aa,"");
        if(verb) sf_raxa(aa);
        sf_oaxa(Fe,aa,3);

        aa=sf_maxa(nht2,-nht*dht,dht);
        sf_setlabel(aa,"ht"); sf_setunit(aa,"s");
        if(verb) sf_raxa(aa);
        sf_oaxa(Fe,aa,4);

        sf_oaxa(Fe,ac,5);
      }
      else{
        Fi = sf_output("out"); sf_settype(Fi,SF_COMPLEX);
        sf_oaxa(Fi,amx,1);
        sf_oaxa(Fi,amy,2);
        sf_oaxa(Fi,az, 3);
        sf_oaxa(Fi,ae, 4);
      }
    }

    /* initialize CIP gathers */
    cip = wexcip_init(cub,nhx,nhy,nhz,nht,nhx2,nhy2,nhz2,nht2,nc,dht,oht,Fc,eic);

    /*------------------------------------------------------------*/
    /* I.C. */
    if(all){
      wexcip_for(cub,cip,Fws,Fwr,Fi,0,0,0);
      if(drv)
        wexcip_for_drv(cub,cip,Fws,Fwr,Fe);
      else
        wexcip_for(cub,cip,Fws,Fwr,Fe,1,0,0);
    }
    else{
      if(eic){
        if(drv){
          sf_warning(" output EIC derivative !");
          wexcip_for_drv(cub,cip,Fws,Fwr,Fe);
        }
        else{
          sf_warning(" output EIC !");
          wexcip_for(cub,cip,Fws,Fwr,Fe,1,0,0);
        }
      }
      else
        wexcip_for(cub,cip,Fws,Fwr,Fi,0,0,0);
    }

    /*------------------------------------------------------------*/
    /* close structures */
    if(all||eic) wexcip_close(cip,1);

    /*------------------------------------------------------------*/
    /* close files */
    if (Fc!=NULL) sf_fileclose(Fc);
    if (Fe!=NULL) sf_fileclose(Fe);
    if (Fi!=NULL) sf_fileclose(Fi);
    if (Fws!=NULL) sf_fileclose(Fws);
    if (Fwr!=NULL) sf_fileclose(Fwr);

    exit (0);
}


