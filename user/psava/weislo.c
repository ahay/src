/* Computing reference slownesses */

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

#include "weiutl.h"
/*^*/

#include "weislo.h"

#define XOOP(a) for(iz=0;  iz<sf_n(cub->az); iz++){ \
                for(ily=0;ily<sf_n(cub->aly);ily++){ \
                for(ilx=0;ilx<sf_n(cub->alx);ilx++){ \
                    {a}                              \
                }}}

/*------------------------------------------------------------*/
weislo3d weislo_init(weicub3d cub,
		     sf_file  Fslo  /* slowness file */
    )
/*< initialize slowness >*/
{
    int iz,jj;
    int nrmax;
    weislo3d slo;

    if (!sf_getint(  "nrmax",&nrmax )) nrmax =     1; /* maximum references */
    /*------------------------------------------------------------*/
    slo = (weislo3d) sf_alloc(1,sizeof(*slo));

    slo->F=Fslo;
    slo->nrmax=nrmax;
    slo->dsmax=cub->dsmax;

    slo->s  = sf_floatalloc3(sf_n(cub->alx),sf_n(cub->aly),sf_n(cub->az));
    slo->sm = sf_floatalloc2(slo->nrmax,sf_n(cub->az));  /* ref slowness squared */
    slo->nr = sf_intalloc              (sf_n(cub->az));  /* nr of ref slownesses @ z */
    
    sf_seek(slo->F,0,SEEK_SET);
    sf_floatread(slo->s[0][0],sf_n(cub->alx)*sf_n(cub->aly)*sf_n(cub->az),slo->F);

    for(iz=0; iz<sf_n(cub->az); iz++)
      for(jj=0; jj<slo->nrmax; jj++)
        slo->sm[iz][jj] = 0.0;

    for (iz=0; iz<sf_n(cub->az); iz++) {
	slo->nr[iz] = weislo(slo->nrmax,
			     slo->dsmax,
			     sf_n(cub->alx)*sf_n(cub->aly),
			     slo->s[iz][0],
			     slo->sm[iz]);
    }
    for (iz=0; iz<sf_n(cub->az)-1; iz++) {
	for (jj=0; jj<slo->nr[iz]; jj++) {
	    slo->sm[iz][jj] = 0.5*(slo->sm[iz][jj]+slo->sm[iz+1][jj]);
	}
    }
    
    return slo;
}

/*------------------------------------------------------------*/
weislo3d weizoslo_init(weicub3d cub,
                     sf_file  Fslo  /* slowness file */
    )
/*< initialize slowness >*/
{
    int iz,jj, ily, ilx;
    int nrmax;
    weislo3d slo;

    if (!sf_getint(  "nrmax",&nrmax )) nrmax =     1; /* maximum references */
    /*------------------------------------------------------------*/
    slo = (weislo3d) sf_alloc(1,sizeof(*slo));

    slo->F=Fslo;
    slo->nrmax=nrmax;
    slo->dsmax=cub->dsmax;
    slo->s  = sf_floatalloc3(sf_n(cub->alx),sf_n(cub->aly),sf_n(cub->az));
    slo->sm = sf_floatalloc2(slo->nrmax,sf_n(cub->az));  /* ref slowness squared */
    slo->nr = sf_intalloc              (sf_n(cub->az));  /* nr of ref slownesses @ z */

    sf_seek(slo->F,0,SEEK_SET);
    sf_floatread(slo->s[0][0],sf_n(cub->alx)*sf_n(cub->aly)*sf_n(cub->az),slo->F);
    XOOP( slo->s[iz][ily][ilx] *= 2.0; );

    for (iz=0; iz<sf_n(cub->az); iz++) {
        slo->nr[iz] = weislo(slo->nrmax,
                             slo->dsmax,
                             sf_n(cub->alx)*sf_n(cub->aly),
                             slo->s[iz][0],
                             slo->sm[iz]);
    }
    for (iz=0; iz<sf_n(cub->az)-1; iz++) {
        for (jj=0; jj<slo->nr[iz]; jj++) {
            slo->sm[iz][jj] = 0.5*(slo->sm[iz][jj]+slo->sm[iz+1][jj]);
        }
    }
    return slo;
}

/*------------------------------------------------------------*/
void weislo_inp(weicub3d cub,
		sf_file Fslo)
/*< input slow dimensions >*/ 
{
    sf_axis alx,aly,az; /* slow axes */

    alx = sf_iaxa(Fslo,1); 
    aly = sf_iaxa(Fslo,2);
    az  = sf_iaxa(Fslo,3);

    sf_copyaxis(cub->alx,alx); sf_setlabel(cub->alx,"lx");
    sf_copyaxis(cub->aly,aly); sf_setlabel(cub->aly,"ly");
    sf_copyaxis(cub->az, az ); sf_setlabel(cub->az,  "z");
}

/*------------------------------------------------------------*/
void weislo_report(weislo3d slo)
/*< report slo parameters >*/ 
{
    sf_warning("nrmax=%d",slo->nrmax);
}

/*------------------------------------------------------------*/
int weislo(int nr           /* maximum number of references */, 
	  float ds         /* minimum slowness separation */, 
	  int ns           /* number of slownesses */, 
	  const float* ss  /* [ns] slowness array */, 
	  float* sr        /* [nr] reference slownesses squared */) 
/*< compute reference slownesses, return their number >*/
{
    int is,jr,ir;
    float smin, smax, s, s2=0., qr, *ss2;

    ss2 = sf_floatalloc(ns);
    for (is=0; is<ns; is++) {
	ss2[is]=ss[is];
    }
    
    smax = sf_quantile(ns-1,ns,ss2);
    smin = sf_quantile(   0,ns,ss2);
    nr = SF_MIN(nr,1+(smax-smin)/ds);
    
    jr=0;
    for (ir=0; ir<nr; ir++) {
	qr = (ir+1.0)/nr - 0.5 * 1./nr;
	s = sf_quantile(qr*ns,ns,ss2);
	if (0==ir || SF_ABS(s-s2) > ds) {
	    sr [jr] = s*s;
	    s2 = s;
	    jr++;
	}
    }
    
    free(ss2);
    return jr;
}

/*------------------------------------------------------------*/
void weislo_close(weislo3d slo)
/*< close slowness >*/
{
    free(**slo->s);  free( *slo->s);  free( slo->s);
    ;                free( *slo->sm); free( slo->sm);
    ;                                 free( slo->nr);
}

