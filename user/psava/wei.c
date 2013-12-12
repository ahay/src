/* 3-D SSR extrapolation using extended split-step */

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

#include <math.h>

#include <rsf.h>
/*^*/

#ifdef _OPENMP
#include <omp.h>
#endif

#include "wei.h"
#include "weitap.h"
#include "weislo.h"
#include "weissr.h"

#include "weiutl.h"
/*^*/

#ifdef SF_HAS_COMPLEX_H
#define cWGH(a,b,c) ((conjf(a)*b)*c)
#define cAWGH(a,b,c) ((a*b)*c)
#else
#define cWGH(a,b,c) sf_cmul(sf_cmul(conjf(a),b),c)
#define cAWGH(a,b,c) sf_cmul(sf_cmul(a,b),c)
#endif

#define cACOR(a,b)   (     (a)*b)
#define  cCOR(a,b)   (conjf(a)*b)

#ifdef SF_HAS_COMPLEX_H
#define cHIC(a,b,w)   crealf( ((conjf(a)*b)*w) )
#define dHIC(a,b,w,e) crealf( ((conjf(a)*b)*w)/(conj(a)*a + e) )
#else
#define cHIC(a,b,w)   crealf(sf_cmul(sf_cmul(conjf(a),b),w))
#define dHIC(a,b,w,e) (crealf(sf_cmul(sf_cmul(conjf(a),b),w))/(sf_crealf(sf_cmul(conjf(a),a)) + e))
#endif

#ifdef SF_HAS_COMPLEX_H
#define cCIC(a,b)   crealf( (conjf(a)*b) )
#define dCIC(a,b,e) crealf( (conjf(a)*b)/(conj(a)*a + e) )
#else
#define cCIC(a,b)   crealf( (sf_cmul(conjf(a),b)) )
#define dCIC(a,b,e) (crealf(sf_cmul(conjf(a),b))/(crealf(sf_cmul(conjf(a),a)) + e))
#endif

#define  YXLOOP(a) \
    for(iy=0;iy<sf_n(cub->amy);iy++){			\
	for(ix=0;ix<sf_n(cub->amx);ix++){		\
	    {a} }} /*  yx loop */

#define CICLOOP(a) for(iz=0;iz<nz;iz++){		\
	for(iy=0;iy<sf_n(cub->amy);iy++){		\
	    for(ix=0;ix<sf_n(cub->amx);ix++){		\
		{a} }}} /* zyx loop */

#define EICLOOP(a) for(ic=0;  ic<nc; ic++ ) {			\
	for(iht=0;iht<sf_n(cub->aht);iht++) {			\
	    for(ihz=0;ihz<sf_n(cub->ahz);ihz++) {		\
		for(ihy=0;ihy<sf_n(cub->ahy);ihy++) {		\
		    for(ihx=0;ihx<sf_n(cub->ahx);ihx++) {	\
			{a} }}}}} /* cc-ht-hz-hy-hx loop */

#define HICLOOP(a) for(ic=0;  ic<nc; ic++ ) {			\
	for(iht=0;iht<sf_n(cub->aht);iht++) {			\
	    for(ihy=0;ihy<sf_n(cub->ahy);ihy++) {		\
		for(ihx=0;ihx<sf_n(cub->ahx);ihx++) {		\
		    {a} }}}} /* cc-ht-hy-hx loop */

/*------------------------------------------------------------*/
weicub3d wei_cube(bool verb_,
		  int   ompnth_)
/*< initialize SR migration space >*/
{
    weicub3d  cub;
    int   pmx,pmy;      /* padding in the k domain */
    int   tmx,tmy;      /* boundary taper size */
    float     eps;      /* dip filter constant */  
    float   dtmax;      /*     time sensitivity */

    cub = (weicub3d) sf_alloc(1,sizeof(*cub));
    
    if (!sf_getfloat(  "eps",&eps   ))   eps =  0.01; /* stability parameter */
    if (!sf_getfloat("dtmax",&dtmax )) dtmax = 0.004; /* max time error */

    cub->verb   = verb_;
    cub->ompnth = ompnth_;

    cub->amx=sf_maxa(1,0,1); sf_setlabel(cub->amx,"mx");
    cub->amy=sf_maxa(1,0,1); sf_setlabel(cub->amy,"my");
    cub->af =sf_maxa(1,0,1); sf_setlabel(cub->af ,"f "); /* frequency in Hertz */
    cub->aw =sf_maxa(1,0,1); sf_setlabel(cub->aw ,"w "); /* frequency in radians */
    cub->ae =sf_maxa(1,0,1); sf_setlabel(cub->ae ,"e "); /* experiments */

    cub->ac =sf_maxa(1,0,1); sf_setlabel(cub->ac, "c ");
    cub->ahx=sf_maxa(1,0,1); sf_setlabel(cub->ahx,"hx");
    cub->ahy=sf_maxa(1,0,1); sf_setlabel(cub->ahy,"hy");
    cub->ahz=sf_maxa(1,0,1); sf_setlabel(cub->ahz,"hz");
    cub->aht=sf_maxa(1,0,1); sf_setlabel(cub->aht,"ht");

    cub->alx=sf_maxa(1,0,1); sf_setlabel(cub->alx,"lx");
    cub->aly=sf_maxa(1,0,1); sf_setlabel(cub->alx,"ly");
    cub->az =sf_maxa(1,0,1); sf_setlabel(cub->az ,"z ");
 
    if (!sf_getint("pmx",&pmx)) pmx=0; /* padding on x */
    if (!sf_getint("pmy",&pmy)) pmy=0; /* padding on y */
    if (!sf_getint("tmx",&tmx)) tmx=0; /*   taper on x */
    if (!sf_getint("tmy",&tmy)) tmy=0; /*   taper on y */
    cub->pmx=pmx;
    cub->pmy=pmy;
    cub->tmx=tmx;
    cub->tmy=tmy;
    
    cub->eps   = eps;
    cub->dtmax = dtmax;
    cub->dsmax = cub->dtmax/sf_d(cub->az);

    return cub;
}

/*------------------------------------------------------------*/
void wei_report(weicub3d cub)
/*< report cube >*/
{
    sf_raxa(cub->amx);
    sf_raxa(cub->amy);
    sf_raxa(cub->af);
    sf_raxa(cub->ae);
    
    sf_raxa(cub->alx);
    sf_raxa(cub->aly);
    sf_raxa(cub->az);

    sf_raxa(cub->ahx);
    sf_raxa(cub->ahy);
    sf_raxa(cub->ahz);
    sf_raxa(cub->aht);
}

/*------------------------------------------------------------*/
void weidat_inp(weicub3d cub,
		sf_file Fdat)
/*< input data dimensions >*/ 
{
    sf_axis amx,amy,af,ae; /* data axes */

    amx = sf_iaxa(Fdat,1); 
    amy = sf_iaxa(Fdat,2);
    af  = sf_iaxa(Fdat,3);
    ae  = sf_iaxa(Fdat,4);

    sf_copyaxis(cub->amx,amx); sf_setlabel(cub->amx,"mx");
    sf_copyaxis(cub->amy,amy); sf_setlabel(cub->amy,"my");
    cub->nxy= sf_n(cub->amx)*sf_n(cub->amy);

    /* read frequency in HZ */
    sf_copyaxis(cub->af, af);  sf_setlabel(cub->af, "f" ); sf_setunit(cub->af,"Hz");
    /* transform frequency to radians */
    sf_setn(cub->aw,         sf_n(af));
    sf_setd(cub->aw,2.*SF_PI*sf_d(af));
    sf_seto(cub->aw,2.*SF_PI*sf_o(af));
    sf_setlabel(cub->aw,"w");
    sf_setunit (cub->aw," ");

    sf_copyaxis(cub->ae, ae ); sf_setlabel(cub->ae , "e");
}

/*------------------------------------------------------------*/
void gensou_inp(weicub3d cub,
                sf_file Fdat)
/*< input generailized source dimensions >*/
{
    sf_axis ac,af; /* data axes */

    ac = sf_iaxa(Fdat,1);
    af = sf_iaxa(Fdat,2);

    sf_copyaxis(cub->amx,cub->alx); sf_setlabel(cub->amx,"mx");
    sf_copyaxis(cub->amy,cub->aly); sf_setlabel(cub->amy,"my");
    sf_copyaxis(cub->af,af); sf_setlabel(cub->af, "f" ); sf_setunit(cub->af,"Hz");
    sf_copyaxis(cub->ac,ac); sf_setlabel(cub->ac, "cc"); sf_setunit(cub->ac,""  );

    /* frequency in radians */
    sf_setn(cub->aw,         sf_n(af));
    sf_setd(cub->aw,2.*SF_PI*sf_d(af));
    sf_seto(cub->aw,2.*SF_PI*sf_o(af));
    sf_setlabel(cub->aw,"w");
    sf_setunit (cub->aw," ");

    cub->nxy= sf_n(cub->amx)*sf_n(cub->amy);
}

/*------------------------------------------------------------*/
void weieic_inp(weicub3d cub,
		sf_file Fcoo)
/*< input EIC parameters >*/ 
{
    sf_axis ac;
    int nhx,nhy,nhz,nht;
    float dht;

    ac  = sf_iaxa(Fcoo,2);

    sf_copyaxis(cub->ac,ac); sf_setlabel(cub->ac,"c"); sf_setunit (cub->ac,"");

    if(! sf_getint  ("nhx",&nhx)) nhx=0; /* x lags */
    if(! sf_getint  ("nhy",&nhy)) nhy=0; /* y lags */
    if(! sf_getint  ("nhz",&nhz)) nhz=0; /* z lags */
    if(! sf_getint  ("nht",&nht)) nht=0; /* t lags */
    if(! sf_getfloat("dht",&dht)) sf_error("need dht");

    cub->ahx=sf_maxa(2*nhx+1,-nhx*sf_d(cub->amx),sf_d(cub->amx)); 
    sf_setlabel(cub->ahx,"hx");

    cub->ahy=sf_maxa(2*nhy+1,-nhy*sf_d(cub->amy),sf_d(cub->amy)); 
    sf_setlabel(cub->ahy,"hy");

    cub->ahz=sf_maxa(2*nhz+1,-nhz*sf_d(cub->az ),sf_d(cub->az )); 
    sf_setlabel(cub->ahz,"hz");

    cub->aht=sf_maxa(2*nht+1,-nht*dht,dht);            
    sf_setlabel(cub->aht,"ht"); 
    sf_setunit (cub->aht,"s");
}

/*------------------------------------------------------------*/
void weihic_inp(weicub3d cub,
		sf_file Fcoo)
/*< input HIC parameters >*/ 
{
    weieic_inp(cub,Fcoo);

    /* no lags on the z axis */
    sf_setn(cub->ahz,1);
    sf_seto(cub->ahz,0.0);
}

/*------------------------------------------------------------*/
void weiwfl_inp(weicub3d cub,
                sf_file Fwfl)
/*< input wavefield dimensions >*/
{
    sf_axis amx,amy,az,af; /* data axes */

    amx = sf_iaxa(Fwfl,1);
    amy = sf_iaxa(Fwfl,2);
    az  = sf_iaxa(Fwfl,3);
    af  = sf_iaxa(Fwfl,4);

    sf_copyaxis(cub->az,az);   sf_setlabel(cub->amy,"z");
    sf_copyaxis(cub->amx,amx); sf_setlabel(cub->amx,"mx");
    sf_copyaxis(cub->amy,amy); sf_setlabel(cub->amy,"my");
    sf_copyaxis(cub->af, af);  sf_setlabel(cub->af, "f" ); sf_setunit(cub->af,"Hz");

    /* frequency in radians */
    sf_setn(cub->aw,         sf_n(af));
    sf_setd(cub->aw,2.*SF_PI*sf_d(af));
    sf_seto(cub->aw,2.*SF_PI*sf_o(af));
    sf_setlabel(cub->aw,"w");
    sf_setunit (cub->aw," ");

    cub->nxy= sf_n(cub->amx)*sf_n(cub->amy);
}


/*------------------------------------------------------------*/
void weiwfl_out(weicub3d cub,
		sf_file Fwfl)
/*< output wfl >*/ 
{
    sf_oaxa(Fwfl,cub->amx,1);
    sf_oaxa(Fwfl,cub->amy,2);
    sf_oaxa(Fwfl,cub->az, 3);
    sf_oaxa(Fwfl,cub->af, 4);
    sf_oaxa(Fwfl,cub->ae, 5);
}

/*------------------------------------------------------------*/
void weidat_out(weicub3d cub,
		sf_file Fwfl)
/*< output dat >*/ 
{
    sf_oaxa(Fwfl,cub->amx,1);
    sf_oaxa(Fwfl,cub->amy,2);
    sf_oaxa(Fwfl,cub->af, 3);
    sf_oaxa(Fwfl,cub->ae, 4);
}

/*------------------------------------------------------------*/
void weicic_out(weicub3d cub,
		sf_file Fcic)
/*< output cic >*/ 
{
    sf_oaxa(Fcic,cub->amx,1);
    sf_oaxa(Fcic,cub->amy,2);
    sf_oaxa(Fcic,cub->az, 3);
}

/*------------------------------------------------------------*/
void weific_out(weicub3d cub,
		sf_file Fcic)
/*< output cic >*/ 
{
    sf_oaxa(Fcic,cub->amx,1);
    sf_oaxa(Fcic,cub->amy,2);
    sf_oaxa(Fcic,cub->az, 3);
    sf_oaxa(Fcic,cub->af, 4);
}


/*------------------------------------------------------------*/
void weieic_out(weicub3d cub,
		sf_file Feic)
/*< output EIC >*/ 
{
    sf_oaxa(Feic,cub->ahx,1);
    sf_oaxa(Feic,cub->ahy,2);
    sf_oaxa(Feic,cub->ahz,3);
    sf_oaxa(Feic,cub->aht,4);
    sf_oaxa(Feic,cub->ac, 5);
}

/*------------------------------------------------------------*/
weiop3d weiwfl_init(weicub3d cub)
/*< initialize weiwfl >*/
{
    weiop3d weop;
    weop = (weiop3d) sf_alloc(1,sizeof(*weop));
    weop->wfld = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    return weop;
}

/*------------------------------------------------------------*/
weiop3d weidtm_init(weicub3d cub)
/*< initialize weidtm >*/
{
    weiop3d weop;
    weop = (weiop3d) sf_alloc(1,sizeof(*weop));
    weop->wfld = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    return weop;
}

/*------------------------------------------------------------*/
weiop3d weicic_init(weicub3d cub)
/*< initialize CIC mig >*/
{
    bool  dflg;
    float deps;

    weiop3d weop;
    weop = (weiop3d) sf_alloc(1,sizeof(*weop));

    /* deconvolution flag and eps */
    if (!sf_getbool( "dflg",&dflg)) dflg = false;
    if (!sf_getfloat("deps",&deps)) deps = 0.1;
    cub->dflg = dflg;
    cub->deps = deps;

    weop->swfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    weop->rwfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    weop->icic = sf_floatalloc3  (sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az));
    return weop;
}

/*------------------------------------------------------------*/
weiop3f weific_init(weicub3d cub)
/*< initialize FIC mig >*/
{
    weiop3f weop;
    weop = (weiop3f) sf_alloc(1,sizeof(*weop));
    weop->swfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    weop->rwfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    weop->icic = sf_floatalloc4  (sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az),sf_n(cub->aw));
    return weop;
}

/*------------------------------------------------------------*/
weiop3d weihic_init(weicub3d cub)
/*< initialize HIC mig >*/
{
    bool  dflg;
    float deps;
    
    weiop3d weop;
    weop = (weiop3d) sf_alloc(1,sizeof(*weop));
    
    /* deconvolution flag and eps */
    if (!sf_getbool( "dflg",&dflg)) dflg = false;
    if (!sf_getfloat("deps",&deps)) deps = 0.1;
    cub->dflg = dflg;
    cub->deps = deps;
    
    weop->swfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    weop->rwfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);

    weop->icic = sf_floatalloc3  (sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az));
    weop->ihic = sf_floatalloc4  (sf_n(cub->ahx),sf_n(cub->ahy),sf_n(cub->aht),sf_n(cub->ac));
    return weop;
}

/*------------------------------------------------------------*/
weiop3d weieic_init(weicub3d cub)
/*< initialize EIC mig >*/
{
    weiop3d weop;
    weop = (weiop3d) sf_alloc(1,sizeof(*weop));

    weop->swfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    weop->rwfl = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);

    weop->icic = sf_floatalloc3  (sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az));
    weop->sone = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az));
    weop->rone = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az));
    weop->ieic = sf_floatalloc5  (sf_n(cub->ahx),sf_n(cub->ahy),sf_n(cub->ahz),sf_n(cub->aht),sf_n(cub->ac));
    return weop;
}

/*------------------------------------------------------------*/
void weiwfl_close(weiop3d weop)
/*< free weiwfl storage >*/
{
    free( **weop->wfld);free(*weop->wfld);free( weop->wfld);
}

/*------------------------------------------------------------*/
void weidtm_close(weiop3d weop)
/*< free weidtm storage >*/
{
    free( **weop->wfld);free(*weop->wfld);free( weop->wfld);
}

/*------------------------------------------------------------*/
void weicic_close(weiop3d weop)
/*< free CIC mig storage >*/
{
    free( **weop->swfl);free(*weop->swfl);free( weop->swfl);
    free( **weop->rwfl);free(*weop->rwfl);free( weop->rwfl);
    free( **weop->icic);free(*weop->icic);free( weop->icic);
}

/*------------------------------------------------------------*/
void weific_close(weiop3f weop)
/*< free FIC mig storage >*/
{
    free( **weop->swfl);free(*weop->swfl);free( weop->swfl);
    free( **weop->rwfl);free(*weop->rwfl);free( weop->rwfl);
    free( ***weop->icic);
    free( **weop->icic);free(*weop->icic);free( weop->icic);
}

/*------------------------------------------------------------*/
void weihic_close(weiop3d weop)
/*< free HIC mig storage >*/
{
    ;                   free(**weop->swfl);free(*weop->swfl);free( weop->swfl);
    ;                   free(**weop->rwfl);free(*weop->rwfl);free( weop->rwfl);
    ;                   free(**weop->icic);free(*weop->icic);free( weop->icic);
    free(***weop->ihic);free(**weop->ihic);free(*weop->ihic);free(weop->ihic);
}

/*------------------------------------------------------------*/
void weieic_close(weiop3d weop)
/*< free EIC mig storage >*/
{
    free(**weop->swfl);free(*weop->swfl);free( weop->swfl);
    free(**weop->rwfl);free(*weop->rwfl);free( weop->rwfl);
    free(**weop->icic);free(*weop->icic);free( weop->icic);

    free(**weop->sone);free(*weop->sone);free( weop->sone);
    free(**weop->rone);free(*weop->rone);free( weop->rone);
    free(****weop->ieic);free(***weop->ieic);free(**weop->ieic);free(*weop->ieic);free(weop->ieic);
}

/*------------------------------------------------------------*/
void weiwfl(weiop3d weop,
	    weicub3d cub,
	    weissr3d ssr,
	    weitap3d tap,
	    weislo3d slo,
	    sf_file Fdat,
	    sf_file Fwfl,
	    bool causal)
/*< wei wavefield >*/
{
    int iz,iw, nw;
    sf_complex w;
    int ompith=0;
    int weisign;

    weisign=causal?+1:-1;

    nw = sf_n(cub->aw);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,w,iz)			\
    shared(Fdat,Fwfl,weop,cub,ssr,tap,slo,nw)
#endif
    for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
	ompith = omp_get_thread_num();
#endif

	/* frequency */
	w = sf_cmplx( cub->eps*sf_d(cub->aw), weisign * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	    if(cub->verb) sf_warning ("WFL ... (ith=%2d) ... <iw=%3d of %3d>",
				      ompith,iw+1,sf_n(cub->aw));
	    
	    /* read data */
	    sf_seek(Fdat,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
	    sf_complexread(weop->wfld[ompith][0],cub->nxy,Fdat);
	    
	    /* write wavefield */
	    iz=0; /* @ the surface */
	    sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
	    sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);	    
	}

	/*------------------------------------------------------------*/
	/* loop over z */
	for (iz=1; iz<sf_n(cub->az); iz++) {

	    weissr(w,weop->wfld[ompith],cub,ssr,slo,iz-1,ompith); /* wei   */
	    weitap(  weop->wfld[ompith],tap);                     /* taper */
	    
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
		/* output wavefield */
		sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
		sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
	    }
	    
	} /* z */
	/*------------------------------------------------------------*/
	
    } /* w */
}

/*------------------------------------------------------------*/
void genwfl(weiop3d weop,
            weicub3d cub,
            weissr3d ssr,
            weitap3d tap,
            weislo3d slo,
            weico3d eico,
            sf_file Fsou,
            sf_file Fwfl,
            bool    down,
            bool   causal)
/*< wei wavefield >*/
{
    int iz,ix,iy,iw,nw,ompith=0,weisign;
    sf_complex w;
    sf_complex **gensou;

    gensou = sf_complexalloc2(sf_n(cub->ac),sf_n(cub->aw));
    weisign=causal?+1:-1;

    /*------------------------------------------------------------*/
    /* read generalized source  */
    sf_complexread(gensou[0],sf_n(cub->ac)*sf_n(cub->aw),Fsou);

    nw = sf_n(cub->aw);

    if(down){ /* downward continuation */

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,w,iz,ix,iy)		\
    shared(Fwfl,cub,ssr,tap,slo,eico,gensou)
#endif
        for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
            ompith = omp_get_thread_num();
#endif
            YXLOOP( weop->wfld[ompith][iy][ix] = sf_cmplx(0.0,0.0); );

            /* frequency */
            w = sf_cmplx( cub->eps*sf_d(cub->aw), weisign * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if(cub->verb) sf_warning ("DC WFL ... (ith=%2d) ... <iw=%3d of %3d>",
					  ompith,iw+1,sf_n(cub->aw));

                /* wavefield injection */
                gensou_inject(weop->wfld[ompith],gensou,eico,0,iw);

                /* write wavefield */
                iz=0; /* @ the surface */
                sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
                sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
            }

            /*------------------------------------------------------------*/
            /* loop over z */
            for (iz=0; iz<sf_n(cub->az)-1; iz++) {
                weissr1(w,weop->wfld[ompith],cub,ssr,slo,iz,ompith,down); /* wei   */
                weitap(   weop->wfld[ompith],tap);                        /* taper */
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    /* wavefield injection */
                    gensou_inject(weop->wfld[ompith],gensou,eico,iz+1,iw);

                    /* output wavefield */
                    sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz+1),SEEK_SET);
                    sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
                }
            } /* z */
        } /* w */
    }
    else{ /* upward continuation */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,w,iz,ix,iy)		\
    shared(Fwfl,cub,ssr,tap,slo,eico,gensou)
#endif
        for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
            ompith = omp_get_thread_num();
#endif
            YXLOOP( weop->wfld[ompith][iy][ix] = sf_cmplx(0.0,0.0); );
            /* frequency */
            w = sf_cmplx( cub->eps*sf_d(cub->aw), weisign * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if(cub->verb) sf_warning ("UC WFL ... (ith=%2d) ... <iw=%3d of %3d>",
					  ompith,iw+1,sf_n(cub->aw));

                iz=sf_n(cub->az)-1; /* @ the bottom */
                /* wavefield injection */
                gensou_inject(weop->wfld[ompith],gensou,eico,iz,iw);

                /* write wavefield */
                sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
                sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
            }
            /*------------------------------------------------------------*/
            /* loop over z */
            for (iz=sf_n(cub->az)-1; iz>0; iz--) {
                weissr1(w,weop->wfld[ompith],cub,ssr,slo,iz,ompith,down);   /* wei   */
                weitap(   weop->wfld[ompith],tap);                          /* taper */
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    /* wavefield injection */
                    gensou_inject(weop->wfld[ompith],gensou,eico,iz-1,iw);

                    /* output wavefield */
                    sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
                    sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
                }
            } /* z */
        } /* w */
    }
}

/*------------------------------------------------------------*/
void adjwfl(weiop3d weop,
            weicub3d cub,
            weissr3d ssr,
            weitap3d tap,
            weislo3d slo,
            sf_file Fsou,
            sf_file Fwfl,
            bool    down,
            bool   causal)
/*< wei wavefield >*/
{
    int iz,ix,iy,iw,nw,ompith=0,weisign;
    sf_complex w;
    sf_complex ****sou;

    sou = sf_complexalloc4(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az),cub->ompnth);
    weisign=causal?+1:-1;

    nw = sf_n(cub->aw);

    if(down){ /* downward continuation */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,w,iz,ix,iy)		\
    shared(Fwfl,Fsou,cub,ssr,tap,slo,sou)
#endif
        for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
            ompith = omp_get_thread_num();
#endif
            YXLOOP( weop->wfld[ompith][iy][ix] = sf_cmplx(0.0,0.0); );

            /* frequency */
            w = sf_cmplx( cub->eps*sf_d(cub->aw), weisign * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if(cub->verb) sf_warning ("DC WFL ... (ith=%2d) ... <iw=%3d of %3d>",
					  ompith,iw+1,sf_n(cub->aw));

                /* read data */
                sf_seek(Fsou,cub->nxy*sf_n(cub->az)*iw*sizeof(sf_complex),SEEK_SET);
                sf_complexread(sou[ompith][0][0],cub->nxy*sf_n(cub->az),Fsou);
            }

            iz=0; /* @ the surface */
            /* wavefield injection */
#ifdef SF_HAS_COMPLEX_H
            YXLOOP( weop->wfld[ompith][iy][ix] += sou[ompith][iz][iy][ix]; );
#else
	    YXLOOP( weop->wfld[ompith][iy][ix] = sf_cadd(weop->wfld[ompith][iy][ix],sou[ompith][iz][iy][ix]); );
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
                sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
            }

            /*------------------------------------------------------------*/
            /* loop over z */
            for (iz=0; iz<sf_n(cub->az)-1; iz++) {
                weissr1(w,weop->wfld[ompith],cub,ssr,slo,iz,ompith,down); /* wei   */
                weitap(   weop->wfld[ompith],tap);                        /* taper */

                /* wavefield injection */
#ifdef SF_HAS_COMPLEX_H
                YXLOOP( weop->wfld[ompith][iy][ix] += sou[ompith][iz+1][iy][ix]; );
#else
		YXLOOP( weop->wfld[ompith][iy][ix] = sf_cadd(weop->wfld[ompith][iy][ix],sou[ompith][iz+1][iy][ix]); );
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    /* output wavefield */
                    sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz+1),SEEK_SET);
                    sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
                }
            } /* z */
        } /* w */
    }
    else{ /* upward continuation */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,w,iz,ix,iy)		\
    shared(Fwfl,cub,ssr,tap,slo,sou)
#endif
        for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
            ompith = omp_get_thread_num();
#endif
            YXLOOP( weop->wfld[ompith][iy][ix] = sf_cmplx(0.0,0.0); );
            /* frequency */
            w = sf_cmplx( cub->eps*sf_d(cub->aw), weisign * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if(cub->verb) sf_warning ("UC WFL ... (ith=%2d) ... <iw=%3d of %3d>",
					  ompith,iw+1,sf_n(cub->aw));
                /* read data */
                sf_seek(Fsou,cub->nxy*sf_n(cub->az)*iw*sizeof(sf_complex),SEEK_SET);
                sf_complexread(sou[ompith][0][0],cub->nxy*sf_n(cub->az),Fsou);
            }

            iz=sf_n(cub->az)-1; /* @ the bottom */
            /* wavefield injection */
#ifdef SF_HAS_COMPLEX_H
            YXLOOP( weop->wfld[ompith][iy][ix] += sou[ompith][iz][iy][ix]; );
#else
	    YXLOOP( weop->wfld[ompith][iy][ix] = sf_cadd(weop->wfld[ompith][iy][ix],sou[ompith][iz][iy][ix]); );
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                /* write wavefield */
                sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
                sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
            }
            /*------------------------------------------------------------*/
            /* loop over z */
            for (iz=sf_n(cub->az)-1; iz>0; iz--) {
#ifdef _OPENMP
#pragma omp critical
#endif
		{
		    /* write wavefield */
		    sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
		    sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
		}

                weissr1(w,weop->wfld[ompith],cub,ssr,slo,iz,ompith,down);   /* wei   */
                weitap(   weop->wfld[ompith],tap);                          /* taper */


                /* wavefield injection */
#ifdef SF_HAS_COMPLEX_H
                YXLOOP( weop->wfld[ompith][iy][ix] += sou[ompith][iz-1][iy][ix]; );
#else
		YXLOOP( weop->wfld[ompith][iy][ix] = sf_cadd(weop->wfld[ompith][iy][ix],sou[ompith][iz-1][iy][ix]); );
#endif
            } /* z */
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
                /* output wavefield */
                sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*(iw)),SEEK_SET);
                sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
	    }
        } /* w */
    }

    free( ***sou); free( **sou); free( *sou); free( sou);
}

/*------------------------------------------------------------*/
void gradient(weicub3d cub,
              sf_file Fsou,
              sf_file Frec,
              sf_file Fgrd,
	      bool conj)
/*< wei gradient >*/
{
    int iz,nz, ix,iy,iw,nw,ompith=0;
    sf_complex ****sou;
    sf_complex ****rec;
    sf_complex  ***cic;

    sou = sf_complexalloc4(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az),cub->ompnth);
    rec = sf_complexalloc4(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az),cub->ompnth);
    cic = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az));

    nz = sf_n(cub->az);
    nw = sf_n(cub->aw);

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero CIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)       \
    private(ix,iy,iz) shared(cub)
#endif
    CICLOOP(cic[iz][iy][ix]=sf_cmplx(0.0,0.0););
    if(cub->verb) fprintf(stderr,"OK\n");

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero CIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)       \
    private(ix,iy,iz) shared(cub)
#endif
    CICLOOP(cic[iz][iy][ix]=sf_cmplx(0.0,0.0););
    if(cub->verb) fprintf(stderr,"OK\n");

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,iz,ix,iy)			\
    shared(Fsou,Frec,cub,conj)
#endif
    for (iw=0; iw<nw; iw++) {

#ifdef _OPENMP
        ompith = omp_get_thread_num();
#endif

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if(cub->verb) sf_warning ("GRD ... (ith=%2d) ... <iw=%3d of %3d>",
				      ompith,iw+1,sf_n(cub->aw));

            /* read data */
            sf_seek(Fsou,cub->nxy*sf_n(cub->az)*iw*sizeof(sf_complex),SEEK_SET);
            sf_complexread(sou[ompith][0][0],cub->nxy*sf_n(cub->az),Fsou);

            sf_seek(Frec,cub->nxy*sf_n(cub->az)*iw*sizeof(sf_complex),SEEK_SET);
            sf_complexread(rec[ompith][0][0],cub->nxy*sf_n(cub->az),Frec);
        }

#ifdef _OPENMP
#pragma omp critical
#endif
        for_cic(cub,sou[ompith],rec[ompith],cic,conj);
    }

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"write gradient...");
    sf_complexwrite(cic[0][0],sf_n(cub->amx)*sf_n(cub->az),Fgrd);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    free( ***sou); free( **sou); free( *sou); free( sou);
    free( ***rec); free( **rec); free( *rec); free( rec);
    free( **cic); free( *cic); free( cic);

}

/*------------------------------------------------------------*/
void weidtm(weiop3d weop,
	    weicub3d cub,
	    weissr3d ssr,
	    weitap3d tap,
	    weislo3d slo,
	    sf_file Fdat,
	    sf_file Fwfl,
	    bool causal)
/*< wei wavefield >*/
{
    int iz,iw,nw,ie;
    sf_complex w;
    int ompith=0;
    int weisign;

    weisign=causal?+1:-1;

    nw = sf_n(cub->aw);

    for(ie=0; ie<sf_n(cub->ae); ie++){
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,w,iz)			\
    shared(Fdat,Fwfl,weop,cub,ssr,tap,slo,ie)
#endif
	for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
	    ompith = omp_get_thread_num();
#endif
	    
	    /* frequency */
	    w = sf_cmplx( cub->eps*sf_d(cub->aw), weisign * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	    
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
		if(cub->verb) sf_warning ("DTM ... (ith=%2d) ... <iw=%3d of %3d> ... <ie=%3d of %3d>",
					  ompith,iw+1,sf_n(cub->aw),ie+1,sf_n(cub->ae));
		
		/* read data */
		sf_seek(Fdat,cub->nxy*(ie*sf_n(cub->aw)+iw)*sizeof(sf_complex),SEEK_SET);
		sf_complexread(weop->wfld[ompith][0],cub->nxy,Fdat);	    
	    }
	    
	    /*------------------------------------------------------------*/
	    /* loop over z */
	    for (iz=1; iz<sf_n(cub->az); iz++) {
		
		weissr(w,weop->wfld[ompith],cub,ssr,slo,iz-1,ompith); /* wei   */
		weitap(  weop->wfld[ompith],tap);                     /* taper */
		
	    } /* z */
	    /*------------------------------------------------------------*/
	    
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
		/* output dtm data */
		sf_seek(Fwfl,sizeof(sf_complex)*cub->nxy*(sf_n(cub->aw)*ie+iw),SEEK_SET);
		sf_complexwrite(weop->wfld[ompith][0],cub->nxy,Fwfl);
	    }
	    
	} /* w */
    } /* e */
}

/*------------------------------------------------------------*/
void weific(weiop3f weop,
	    weicub3d cub,
	    weissr3d ssr,
	    weitap3d tap,
	    weislo3d slo,
	    sf_file Fsou,
	    sf_file Frec,
	    sf_file Fcic)
/*< FIC migration >*/
{
    int ix,iy,iz,iw,nw;
    sf_complex ws,wr;
    int ompith=0;
   
    nw =  sf_n(cub->aw);

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero CIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(iw,ix,iy,iz) shared(weop,cub,nw)
#endif
    for(iw = 0; iw < nw; ++iw){
	for(iz = 0; iz < sf_n(cub->az); ++iz){
	    for(iy = 0; iy < sf_n(cub->amy); ++iy){
		for(ix = 0; ix < sf_n(cub->amx); ++ix){
		    weop->icic[iw][iz][iy][ix] = 0;
		}
	    }
	}
    }
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,ws,wr,iz)			\
    shared(Fsou,Frec,weop,cub,ssr,tap,slo,nw)
#endif
    for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
	ompith = omp_get_thread_num();
#endif

	/* frequency */
	ws = sf_cmplx( cub->eps*sf_d(cub->aw), (+1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	wr = sf_cmplx( cub->eps*sf_d(cub->aw), (-1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	    if(cub->verb) sf_warning ("WEI ... (ith=%2d) ... <iw=%3d of %3d>",
				      ompith,iw+1,sf_n(cub->aw));
	    
	    /* read data */
	    sf_seek(Fsou,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
	    sf_complexread(weop->swfl[ompith][0],cub->nxy,Fsou);

	    sf_seek(Frec,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
	    sf_complexread(weop->rwfl[ompith][0],cub->nxy,Frec);
	    
	}
	iz=0; /* I.C. @ iz=0 */
	weific_one(weop,cub,ompith,iz,iw);

	/*------------------------------------------------------------*/
	/* loop over z */
	for (iz=1; iz<sf_n(cub->az); iz++) {

	    weissr(ws,weop->swfl[ompith],cub,ssr,slo,iz-1,ompith);
	    weitap(   weop->swfl[ompith],tap);

	    weissr(wr,weop->rwfl[ompith],cub,ssr,slo,iz-1,ompith);
	    weitap(   weop->rwfl[ompith],tap);
	    
	    weific_one(weop,cub,ompith,iz,iw);
	} /* z */
	/*------------------------------------------------------------*/

    } /* w */

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"write CIC image...");
    weific_write(weop,cub,Fcic);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/
}

/*------------------------------------------------------------*/
void weicic(weiop3d weop,
	    weicub3d cub,
	    weissr3d ssr,
	    weitap3d tap,
	    weislo3d slo,
	    sf_file Fsou,
	    sf_file Frec,
	    sf_file Fcic)
/*< CIC migration >*/
{
    int ix,iy,iz,nz,iw,nw;
    sf_complex ws,wr;
    int ompith=0;
   
    nz = sf_n(cub->az);
    nw = sf_n(cub->aw);

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero CIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ix,iy,iz) shared(weop,cub)
#endif
    CICLOOP(weop->icic[iz][iy][ix]=0;);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)	\
    private(ompith,iw,ws,wr,iz)			\
    shared(Fsou,Frec,weop,cub,ssr,tap,slo,nw)
#endif
    for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
	ompith = omp_get_thread_num();
#endif

	/* frequency */
	ws = sf_cmplx( cub->eps*sf_d(cub->aw), (+1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	wr = sf_cmplx( cub->eps*sf_d(cub->aw), (-1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	    if(cub->verb) sf_warning ("WEI ... (ith=%2d) ... <iw=%3d of %3d>",
				      ompith,iw+1,sf_n(cub->aw));
	    
	    /* read data */
	    sf_seek(Fsou,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
	    sf_complexread(weop->swfl[ompith][0],cub->nxy,Fsou);

	    sf_seek(Frec,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
	    sf_complexread(weop->rwfl[ompith][0],cub->nxy,Frec);
	    
	}
	iz=0; /* I.C. @ iz=0 */
	weicic_one(weop,cub,ompith,iz);

	/*------------------------------------------------------------*/
	/* loop over z */
	for (iz=1; iz<sf_n(cub->az); iz++) {

	    weissr(ws,weop->swfl[ompith],cub,ssr,slo,iz-1,ompith);
	    weitap(   weop->swfl[ompith],tap);

	    weissr(wr,weop->rwfl[ompith],cub,ssr,slo,iz-1,ompith);
	    weitap(   weop->rwfl[ompith],tap);
	    
	    weicic_one(weop,cub,ompith,iz);
	} /* z */
	/*------------------------------------------------------------*/

    } /* w */

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"write CIC image...");
    weicic_write(weop,cub,Fcic);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/
}

/*------------------------------------------------------------*/
void weizocic(weiop3d weop,
            weicub3d cub,
            weissr3d ssr,
            weitap3d tap,
            weislo3d slo,
            sf_file Fdat,
            sf_file Fcic)
/*< zero-offset migration >*/
{
    int ix,iy,iz,nz,iw,nw;
    sf_complex ws;
    int ompith=0;

    nz = sf_n(cub->az);
    nw = sf_n(cub->aw);

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero CIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)       \
    private(ix,iy,iz) shared(weop,cub)
#endif
    CICLOOP(weop->icic[iz][iy][ix]=0;);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  \
    private(ompith,iw,ws,iz)                 \
    shared(Fdat,weop,cub,ssr,tap,slo)
#endif
    for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
        ompith = omp_get_thread_num();
#endif

        /* frequency */
        ws = sf_cmplx( cub->eps*sf_d(cub->aw), (-1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if(cub->verb) sf_warning ("WEI ... (ith=%2d) ... <iw=%3d of %3d>",
                                      ompith,iw+1,sf_n(cub->aw));
            /* read data */
            sf_seek(Fdat,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
            sf_complexread(weop->swfl[ompith][0],cub->nxy,Fdat);
        }
        iz=0; /* I.C. @ iz=0 */
        weizocic_one(weop,cub,ompith,iz);

        /*------------------------------------------------------------*/
        /* loop over z */
        for (iz=1; iz<sf_n(cub->az); iz++) {

            weissr(ws,weop->swfl[ompith],cub,ssr,slo,iz-1,ompith);
            weitap(   weop->swfl[ompith],tap);

            weizocic_one(weop,cub,ompith,iz);
        } /* z */
        /*------------------------------------------------------------*/
    } /* w */

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"write CIC image...");
    weicic_write(weop,cub,Fcic);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/
}

/*------------------------------------------------------------*/
void weizomod(weiop3d weop,
            weicub3d cub,
            weissr3d ssr,
            weitap3d tap,
            weislo3d slo,
            sf_file Fdat,
            sf_file Fcic)
/*< zero-offset modeling >*/
{
    int ix,iy,iz,iw,nw;
    sf_complex ws;
    int ompith=0;

    sf_complex ***img;
    img=  sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az));
    sf_complexread(img[0][0],cub->nxy*sf_n(cub->az),Fcic);
    /*------------------------------------------------------------*/

    nw = sf_n(cub->aw);

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)  \
    private(ompith,iw,ws,iz,ix,iy)                   \
    shared(Fdat,weop,cub,ssr,tap,slo)
#endif
    for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
        ompith = omp_get_thread_num();
#endif
        YXLOOP( weop->swfl[ompith][iy][ix] = sf_cmplx(0.0,0.0); );

        /* frequency */
        ws = sf_cmplx( cub->eps*sf_d(cub->aw), (+1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if(cub->verb) sf_warning ("WEI ... (ith=%2d) ... <iw=%3d of %3d>",
                                      ompith,iw+1,sf_n(cub->aw));
        }
        /*------------------------------------------------------------*/
        /* loop over z */
        for (iz=sf_n(cub->az)-1; iz>0; iz--) {
#ifdef SF_HAS_COMPLEX_H
            YXLOOP( weop->swfl[ompith][iy][ix] += img[iz][iy][ix]; );
#else
            YXLOOP( weop->swfl[ompith][iy][ix] = sf_cadd(weop->swfl[ompith][iy][ix],img[iz][iy][ix]); );
#endif
            weitap(   weop->swfl[ompith],tap);
            weissr(ws,weop->swfl[ompith],cub,ssr,slo,iz-1,ompith);
        } /* z */

        /*------------------------------------------------------------*/
#ifdef SF_HAS_COMPLEX_H
        YXLOOP( weop->swfl[ompith][iy][ix] += img[0][iy][ix]; );
#else
        YXLOOP( weop->swfl[ompith][iy][ix] = sf_cadd(weop->swfl[ompith][iy][ix],img[0][iy][ix]); );
#endif
        weitap(   weop->swfl[ompith],tap);
#ifdef _OPENMP
#pragma omp critical
#endif
        {
          sf_seek(Fdat,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
          sf_complexwrite(weop->swfl[ompith][0],cub->nxy,Fdat);
        }
    } /* w */

    /*------------------------------------------------------------*/
    free( **img); free( *img); free( img);
}

/*------------------------------------------------------------*/
weico3d weicoo_init(weicub3d cub,
		    sf_file Fcoo)
/*< EIC coordinates >*/
{
    int ic,icx,icy,icz,ihx,ihy,ihz;
    int nhx,nhy,nhz;
    weico3d eico;
    eico = (weico3d) sf_alloc(1,sizeof(*eico));

    /*------------------------------------------------------------*/
    eico->ccin= sf_boolalloc    (sf_n(cub->ac));

    eico->cc  = (pt3d*) sf_alloc(sf_n(cub->ac),sizeof(*eico->cc));
    pt3dread1(Fcoo,eico->cc,sf_n(cub->ac),3);
    /*------------------------------------------------------------*/
    eico->mcxall=sf_intalloc2(sf_n(cub->ahx),sf_n(cub->ac));
    eico->pcxall=sf_intalloc2(sf_n(cub->ahx),sf_n(cub->ac));

    eico->mcyall=sf_intalloc2(sf_n(cub->ahy),sf_n(cub->ac));
    eico->pcyall=sf_intalloc2(sf_n(cub->ahy),sf_n(cub->ac));

    eico->mczall=sf_intalloc2(sf_n(cub->ahz),sf_n(cub->ac));
    eico->pczall=sf_intalloc2(sf_n(cub->ahz),sf_n(cub->ac));

    eico->iczall=sf_intalloc (sf_n(cub->ac));
    /*------------------------------------------------------------*/
    nhx=(sf_n(cub->ahx)-1)/2;
    nhy=(sf_n(cub->ahy)-1)/2;
    nhz=(sf_n(cub->ahz)-1)/2;

    eico->cxmin = sf_o(cub->amx) +                   nhx *sf_d(cub->amx);
    eico->cxmax = sf_o(cub->amx) + (sf_n(cub->amx)-1-nhx)*sf_d(cub->amx);
    eico->cymin = sf_o(cub->amy) +                   nhy *sf_d(cub->amy);
    eico->cymax = sf_o(cub->amy) + (sf_n(cub->amy)-1-nhy)*sf_d(cub->amy);
    eico->czmin = sf_o(cub->az)  +                   nhz *sf_d(cub->az);
    eico->czmax = sf_o(cub->az)  + (sf_n(cub->az) -1-nhz)*sf_d(cub->az);
    /*------------------------------------------------------------*/
    for(ic=0; ic<sf_n(cub->ac); ic++) {            
	eico->ccin[ic]=(bool) (eico->cc[ic].x>=eico->cxmin && 
			       eico->cc[ic].x<=eico->cxmax &&
			       eico->cc[ic].y>=eico->cymin && 
			       eico->cc[ic].y<=eico->cymax &&
			       eico->cc[ic].z>=eico->czmin && 
			       eico->cc[ic].z<=eico->czmax);
	
	if(eico->ccin[ic]) {
	    
	    icx = 0.5+(eico->cc[ic].x-sf_o(cub->amx))/sf_d(cub->amx);
	    for(ihx=-nhx; ihx<nhx+1; ihx++) {
		eico->mcxall[ic][nhx+ihx] = icx-ihx;
		eico->pcxall[ic][nhx+ihx] = icx+ihx;
	    }
	    
	    icy = 0.5+(eico->cc[ic].y-sf_o(cub->amy))/sf_d(cub->amy);
	    for(ihy=-nhy; ihy<nhy+1; ihy++) {
		eico->mcyall[ic][nhy+ihy] = icy-ihy;
		eico->pcyall[ic][nhy+ihy] = icy+ihy;
	    }
	    
	    icz = 0.5+(eico->cc[ic].z-sf_o(cub->az))/sf_d(cub->az);
	    for(ihz=-nhz; ihz<nhz+1; ihz++) {
		eico->mczall[ic][nhz+ihz] = icz-ihz;
		eico->pczall[ic][nhz+ihz] = icz+ihz;
	    }
	    eico->iczall[ic] = icz;

        } /* end if ccin */

    } /* end cc loop */
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    wei_hicindex(cub,eico);
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    wei_timelag(cub,eico);
    /*------------------------------------------------------------*/

    return eico;
}

/*------------------------------------------------------------*/
weico3d gencoo_init(weicub3d cub,
                    sf_file Fcoo)
/*< EIC coordinates >*/
{
    int ic,icx,icy,icz,nhx,nhy,nhz;
    weico3d eico;
    eico = (weico3d) sf_alloc(1,sizeof(*eico));

    /*------------------------------------------------------------*/
    eico->ccin= sf_boolalloc    (sf_n(cub->ac));
    eico->cc  = (pt3d*) sf_alloc(sf_n(cub->ac),sizeof(*eico->cc));
    pt3dread1(Fcoo,eico->cc,sf_n(cub->ac),3);
    /*------------------------------------------------------------*/
    eico->mcxall=sf_intalloc2(1,sf_n(cub->ac));
    eico->mcyall=sf_intalloc2(1,sf_n(cub->ac));
    eico->mczall=sf_intalloc2(1,sf_n(cub->ac));

    eico->iczall=sf_intalloc (sf_n(cub->ac));
    /*------------------------------------------------------------*/
    nhx = 0;  nhy = 0;  nhz = 0;

    eico->cxmin = sf_o(cub->amx) +                   nhx *sf_d(cub->amx);
    eico->cxmax = sf_o(cub->amx) + (sf_n(cub->amx)-1-nhx)*sf_d(cub->amx);
    eico->cymin = sf_o(cub->amy) +                   nhy *sf_d(cub->amy);
    eico->cymax = sf_o(cub->amy) + (sf_n(cub->amy)-1-nhy)*sf_d(cub->amy);
    eico->czmin = sf_o(cub->az)  +                   nhz *sf_d(cub->az);
    eico->czmax = sf_o(cub->az)  + (sf_n(cub->az) -1-nhz)*sf_d(cub->az);

    /*------------------------------------------------------------*/
    for(ic=0; ic<sf_n(cub->ac); ic++) {
        eico->ccin[ic]=(bool) (eico->cc[ic].x>=eico->cxmin && 
			       eico->cc[ic].x<=eico->cxmax &&
			       eico->cc[ic].y>=eico->cymin && 
			       eico->cc[ic].y<=eico->cymax &&
			       eico->cc[ic].z>=eico->czmin && 
			       eico->cc[ic].z<=eico->czmax);

        if(eico->ccin[ic]) {
            icx = 0.5+(eico->cc[ic].x-sf_o(cub->amx))/sf_d(cub->amx);
            eico->mcxall[ic][0] = icx;

            icy = 0.5+(eico->cc[ic].y-sf_o(cub->amy))/sf_d(cub->amy);
            eico->mcyall[ic][0] = icy;

            icz = 0.5+(eico->cc[ic].z-sf_o(cub->az))/sf_d(cub->az);
            eico->mczall[ic][0] = icz;
            eico->iczall[ic] = icz;
        } /* end if ccin */
    } /* end cc loop */

    /*------------------------------------------------------------*/
    wei_hicindex(cub,eico);

    return eico;
}


/*------------------------------------------------------------*/
void wei_hicindex(weicub3d cub,
		  weico3d eico)
/*< set cc indices >*/
{
    int iz,ic;

    eico->ncofz =sf_intalloc (sf_n(cub->az));
    eico->icofz =(int**) sf_alloc    (sf_n(cub->az),sizeof(int*));

    for(iz=0; iz<sf_n(cub->az); iz++) {

	eico->ncofz[iz]=0;
	for(ic=0; ic<sf_n(cub->ac); ic++) {
	    if(eico->ccin[ic] && eico->iczall[ic] == iz) {
		eico->ncofz[iz]++;
	    }
	}

	eico->icofz[iz] = sf_intalloc(eico->ncofz[iz]+1);	

	eico->ncofz[iz]=0;
	for(ic=0; ic<sf_n(cub->ac); ic++) {
	    if(eico->ccin[ic] && eico->iczall[ic] == iz) {
		eico->icofz[iz][eico->ncofz[iz]] = ic;
		eico->ncofz[iz]++;
	    }
	}

    }
}

/*------------------------------------------------------------*/
void wei_timelag(weicub3d cub,
		 weico3d eico)
/*< time delay init >*/
{
    int iw,iht;
    float ht, ww;

    eico->tt=sf_complexalloc2(sf_n(cub->aht),sf_n(cub->aw));

    for(iw=0; iw<sf_n(cub->aw); iw++) {
	ww = -( sf_o(cub->aw)+iw*sf_d(cub->aw));
	
	for(iht=0; iht<sf_n(cub->aht); iht++) {
	    ht = sf_o(cub->aht) + (iht)*sf_d(cub->aht);
	    eico->tt[iw][iht] = sf_cmplx(cosf(2*ww*ht),sinf(2*ww*ht));
	}
    }

}

/*------------------------------------------------------------*/
void weicoo_close(weico3d eico)
/*< EIC coordinates close >*/
{
    free(*eico->mcxall);free( eico->mcxall);
    free(*eico->pcxall);free( eico->pcxall);

    free(*eico->mcyall);free( eico->mcyall);
    free(*eico->pcyall);free( eico->pcyall);

    free(*eico->mczall);free( eico->mczall);
    free(*eico->pczall);free( eico->pczall);
    
    free( eico->ccin);
    free( eico->cc);

    free(*eico->tt);free( eico->tt);

    free( eico->iczall);
    free(*eico->icofz);free( eico->icofz);
    ;                  free( eico->ncofz);
}

/*------------------------------------------------------------*/
void gencoo_close(weico3d eico)
/*< EIC coordinates close >*/
{
    free(*eico->mcxall);free( eico->mcxall);
    free(*eico->mcyall);free( eico->mcyall);
    free(*eico->mczall);free( eico->mczall);

    free( eico->ccin);
    free( eico->cc);

    free( eico->iczall);
    free(*eico->icofz);free( eico->icofz);
    ;                  free( eico->ncofz);
}

/*------------------------------------------------------------*/
void weihic(weiop3d weop,
	    weicub3d cub,
	    weissr3d ssr,
	    weitap3d tap,
	    weislo3d slo,
	    weico3d eico,
	    sf_file Fsou,
	    sf_file Frec,
	    sf_file Fcic,
	    sf_file Feic)
/*< HIC migration >*/
{
    int ix,iy,iz,nz,iw,nw;
    sf_complex ws,wr;
    int ompith=0;
    int ihx,ihy,iht,ic,nc;
   
    nz = sf_n(cub->az);
    nc = sf_n(cub->ac);
    nw = sf_n(cub->aw);

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero CIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ix,iy,iz) shared(weop,cub)
#endif
    CICLOOP(weop->icic[iz][iy][ix]=0;);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero HIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ihx,ihy,iht,ic) shared(weop,cub)
#endif
    HICLOOP(weop->ihic[ic][iht][ihy][ihx]=0;);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(ompith,iw,ws,wr,iz)			\
    shared(Fsou,Frec,weop,cub,ssr,tap,slo)
#endif
    for (iw=0; iw<nw; iw++) {
#ifdef _OPENMP
	ompith = omp_get_thread_num();
#endif

	/* frequency */
	ws = sf_cmplx( cub->eps*sf_d(cub->aw), (+1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	wr = sf_cmplx( cub->eps*sf_d(cub->aw), (-1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	    if(cub->verb) sf_warning ("WEI ... (ith=%2d) ... <iw=%3d of %3d>",
				      ompith,iw+1,sf_n(cub->aw));
	    
	    /* read data */
	    sf_seek(Fsou,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
	    sf_complexread(weop->swfl[ompith][0],cub->nxy,Fsou);

	    sf_seek(Frec,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
	    sf_complexread(weop->rwfl[ompith][0],cub->nxy,Frec);
	    
	}
	iz=0; /* I.C. @ iz=0 */
	weicic_one(weop,cub,     ompith,iz);
	weihic_one(weop,cub,eico,ompith,iz,iw);

	/*------------------------------------------------------------*/
	/* loop over z */
	for (iz=1; iz<sf_n(cub->az); iz++) {

	    weissr(ws,weop->swfl[ompith],cub,ssr,slo,iz-1,ompith);
	    weitap(   weop->swfl[ompith],tap);

	    weissr(wr,weop->rwfl[ompith],cub,ssr,slo,iz-1,ompith);
	    weitap(   weop->rwfl[ompith],tap);
	    
	    /* I.C. @ iz */
	    weicic_one(weop,cub,     ompith,iz);
	    weihic_one(weop,cub,eico,ompith,iz,iw);
	} /* z */
	/*------------------------------------------------------------*/
    } /* w */

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"write CIC image...");
    weicic_write(weop,cub,Fcic);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"write HIC image...");
    weihic_write(weop,cub,Feic);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/
}

/*------------------------------------------------------------*/
void weieic(weiop3d weop,
	    weicub3d cub,
	    weissr3d ssr,
	    weitap3d tap,
	    weislo3d slo,
	    weico3d eico,
	    sf_file Fsou,
	    sf_file Frec,
	    sf_file Fcic,
	    sf_file Feic)
/*< EIC migration >*/
{
    int ix,iy,iz,nz,iw,ic,nc,ihx,ihy,ihz,iht;
    sf_complex ws,wr;
    int ompith=0;
    int ibk,nbk;
    int jw,nw;
    sf_fslice *stmp,*rtmp;

    nz = sf_n(cub->az);
    nc = sf_n(cub->ac);
 
    /*------------------------------------------------------------*/
    nw=cub->ompnth;       /* frequencies in a block */
    nbk=sf_n(cub->aw)/nw; /* frequency blocks */

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"init slices...");
    stmp = (sf_fslice*) sf_alloc(nw,sizeof(sf_fslice));
    rtmp = (sf_fslice*) sf_alloc(nw,sizeof(sf_fslice));
    for(jw=0;jw<nw;jw++) {
	stmp[jw] = sf_fslice_init(cub->nxy,sf_n(cub->az),sizeof(sf_complex));
	rtmp[jw] = sf_fslice_init(cub->nxy,sf_n(cub->az),sizeof(sf_complex));

	for (iz=0; iz<sf_n(cub->az); iz++) {
	    sf_fslice_put(stmp[jw],iz,weop->swfl[0][0]);
	    sf_fslice_put(rtmp[jw],iz,weop->rwfl[0][0]);
	}
    }
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero CIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ix,iy,iz) shared(weop,cub)
#endif
    CICLOOP(weop->icic[iz][iy][ix]=0;);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"zero EIC image...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)		\
    private(ihx,ihy,ihz,iht,ic) shared(weop,cub)
#endif
    EICLOOP(weop->ieic[ic][iht][ihz][ihy][ihx]=0;);
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    /* loop over blocks */
    for (ibk=0;ibk<nbk;ibk++) {   

	/*------------------------------------------------------------*/ 
	/* W.R. loop */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)				\
    private(ompith,iw,jw,ws,wr,iz)					\
    shared(Fsou,Frec,Feic,Fcic,weop,cub,ssr,tap,slo,eico,nw,ibk)
#endif
	for (jw=0;jw<nw;jw++) {
#ifdef _OPENMP
	    ompith = omp_get_thread_num();
#endif
	    /* frequency */
	    iw=ibk*nw+jw;
	    ws = sf_cmplx( cub->eps*sf_d(cub->aw), (+1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	    wr = sf_cmplx( cub->eps*sf_d(cub->aw), (-1) * (sf_o(cub->aw)+iw*sf_d(cub->aw)) );
	    
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
		if(cub->verb) sf_warning (" WR ... (ith=%2d) ... <iw=%3d of %3d>",
					  ompith,iw+1,sf_n(cub->aw));

		/* read data */
		sf_seek(Fsou,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
		sf_complexread(weop->swfl[ompith][0],cub->nxy,Fsou);
		
		sf_seek(Frec,cub->nxy*iw*sizeof(sf_complex),SEEK_SET);
		sf_complexread(weop->rwfl[ompith][0],cub->nxy,Frec);
	    }
	    /* write wavefield */
	    iz=0;
	    sf_fslice_put(stmp[jw],iz,weop->swfl[ompith][0]);
	    sf_fslice_put(rtmp[jw],iz,weop->rwfl[ompith][0]);

	    /*------------------------------------------------------------*/
	    /* loop over z */
	    for (iz=1; iz<sf_n(cub->az); iz++) {
		
		weissr(ws,weop->swfl[ompith],cub,ssr,slo,iz-1,ompith);
		weitap(   weop->swfl[ompith],tap);

		weissr(wr,weop->rwfl[ompith],cub,ssr,slo,iz-1,ompith);
		weitap(   weop->rwfl[ompith],tap);

		sf_fslice_put(stmp[jw],iz,weop->swfl[ompith][0]);
		sf_fslice_put(rtmp[jw],iz,weop->rwfl[ompith][0]);
	    } /* z */
	    /*------------------------------------------------------------*/
	   	    
	} /* w */
	/*------------------------------------------------------------*/
	
	if(sf_n(cub->ac)>cub->ompnth)
	    sf_warning("IC parallelized on CIP");
	else
	    sf_warning("IC parallelized on tau");

	/*------------------------------------------------------------*/
	/* I.C. loop */
	/*------------------------------------------------------------*/
	for(jw=0;jw<nw;jw++) {
	    sf_warning(" IC ... (jw=%2d)",jw);
	    
	    for (iz=0; iz<sf_n(cub->az); iz++) {
		sf_fslice_get(stmp[jw],iz,weop->sone[iz][0]);
		sf_fslice_get(rtmp[jw],iz,weop->rone[iz][0]);
	    }

	    weicic_apply(weop,cub);
	    weieic_apply(jw,ibk,weop,eico,cub);

	} /* jw */
	/*------------------------------------------------------------*/

	/*------------------------------------------------------------*/
	if(cub->verb) fprintf(stderr,"write CIC image...");
	weicic_write(weop,cub,Fcic);
	if(cub->verb) fprintf(stderr,"OK\n");
	/*------------------------------------------------------------*/
	
	/*------------------------------------------------------------*/
	if(cub->verb) fprintf(stderr,"write EIC image...");
	weieic_write(weop,cub,Feic);
	if(cub->verb) fprintf(stderr,"OK\n");
	/*------------------------------------------------------------*/
	
    } /* nbk */
    
    /*------------------------------------------------------------*/
    if(cub->verb) fprintf(stderr,"close slices...");
    for(jw=0;jw<nw;jw++) {
	sf_fslice_close(stmp[jw]);
	sf_fslice_close(rtmp[jw]);
    }
    if(cub->verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

}


/*------------------------------------------------------------*/
void weicic_write(weiop3d weop,
		  weicub3d cub,
		  sf_file Fcic)
/*< write CIC >*/
{
    sf_seek(Fcic,0,SEEK_SET);
    sf_floatwrite(weop->icic[0][0],
		  cub->nxy*sf_n(cub->az),
		  Fcic);
}

/*------------------------------------------------------------*/
void weific_write(weiop3f weop,
		  weicub3d cub,
		  sf_file Fcic)
/*< write CIC >*/
{
    sf_seek(Fcic,0,SEEK_SET);
    sf_floatwrite(weop->icic[0][0][0],
		  cub->nxy*sf_n(cub->az)*sf_n(cub->aw),
		  Fcic);
}

/*------------------------------------------------------------*/
void weihic_write(weiop3d weop,
		  weicub3d cub,
		  sf_file Feic)
/*< write HIC >*/
{
    sf_seek(Feic,0,SEEK_SET);
    sf_floatwrite(weop->ihic[0][0][0],
		  sf_n(cub->ac)*sf_n(cub->aht)*sf_n(cub->ahy)*sf_n(cub->ahx),
		  Feic);
}


/*------------------------------------------------------------*/
void weieic_write(weiop3d weop,
		  weicub3d cub,
		  sf_file Feic)
/*< write EIC >*/
{
    sf_seek(Feic,0,SEEK_SET);
    sf_floatwrite(weop->ieic[0][0][0][0],
		  sf_n(cub->ac)*sf_n(cub->aht)*sf_n(cub->ahz)*sf_n(cub->ahy)*sf_n(cub->ahx),
		  Feic);
}


/*------------------------------------------------------------*/
void weicic_one(weiop3d weop,
		weicub3d cub,
		int ompith,
		int iz)
/*< apply CIC at one z >*/
{
    int ix,iy;
    float meps;
    
    if(cub->dflg) { /* deconvolution */

/* 
	meps=0;
	for(iy=0;iy<sf_n(cub->amy);iy++){
	    for(ix=0;ix<sf_n(cub->amx);ix++){
		meps = SF_MAX(meps,
			      cCIC(weop->swfl[ompith][iy][ix],
				   weop->swfl[ompith][iy][ix]));
	    }
	}
	meps *= cub->deps;
	meps = SF_MAX(meps,1e-5);
*/
	meps = cub->deps;

	for(iy=0;iy<sf_n(cub->amy);iy++){
	    for(ix=0;ix<sf_n(cub->amx);ix++){
#ifdef _OPENMP
#pragma omp atomic
#endif	
		weop->icic[iz][iy][ix] +=
		    dCIC(weop->swfl[ompith][iy][ix],
			 weop->rwfl[ompith][iy][ix],
			 meps);

	    }
	}
	

	for(iy=0;iy<sf_n(cub->amy);iy++){
	    for(ix=0;ix<sf_n(cub->amx);ix++){
#ifdef _OPENMP
#pragma omp atomic
#endif	
		weop->icic[iz][iy][ix] +=
		    cCIC(weop->swfl[ompith][iy][ix],
			 weop->rwfl[ompith][iy][ix]);
	    }
	}

    } else {        /* correlation */

	for(iy=0;iy<sf_n(cub->amy);iy++){
	    for(ix=0;ix<sf_n(cub->amx);ix++){
#ifdef _OPENMP
#pragma omp atomic
#endif	
		weop->icic[iz][iy][ix] +=
		    cCIC(weop->swfl[ompith][iy][ix],
			 weop->rwfl[ompith][iy][ix]);
	    }
	}
    }
    
}

/*------------------------------------------------------------*/
void weizocic_one(weiop3d weop,
                weicub3d cub,
                int ompith,
                int iz)
/*< apply CIC at one z >*/
{
    int ix,iy;

    for(iy=0;iy<sf_n(cub->amy);iy++)
        for(ix=0;ix<sf_n(cub->amx);ix++)
            weop->icic[iz][iy][ix] += crealf(weop->swfl[ompith][iy][ix]);
}

/*------------------------------------------------------------*/
void weific_one(weiop3f weop,
		weicub3d cub,
		int ompith,
		int iz,
		int iw)
/*< apply FIC at one z >*/
{
    int ix,iy;
    
    for(iy=0;iy<sf_n(cub->amy);iy++){
	for(ix=0;ix<sf_n(cub->amx);ix++){
	    
#ifdef _OPENMP
#pragma omp atomic
#endif	
	    weop->icic[iw][iz][iy][ix] +=
		cCIC(weop->swfl[ompith][iy][ix],
		     weop->rwfl[ompith][iy][ix]);
	    
	}
    }
}

/*------------------------------------------------------------*/
void weihic_one(weiop3d weop,
		weicub3d cub,
		weico3d eico,
		int ompith,
		int iz,
		int iw)
/*< apply HIC at one z >*/
{
    int ihx,ihy,iht,ic,jc;
    sf_complex wt;
    int mcx,pcx,ix;
    int mcy,pcy,iy;
    float meps;

    if(cub->dflg) { /* deconvolution */

	meps=0;
	for(iy=0;iy<sf_n(cub->amy);iy++){
	    for(ix=0;ix<sf_n(cub->amx);ix++){
		meps = SF_MAX(meps,
			      cCIC(weop->swfl[ompith][iy][ix],
				   weop->swfl[ompith][iy][ix]));
	    }
	}
	meps *= cub->deps;
	meps = SF_MAX(meps,1e-5);

	for(jc=0;jc<eico->ncofz[iz];jc++) {
	    ic=eico->icofz[iz][jc];
	    
	    for        (iht=0; iht<sf_n(cub->aht); iht++) { wt =eico->tt[iw][iht];
		for    (ihy=0; ihy<sf_n(cub->ahy); ihy++) { mcy=eico->mcyall[ic][ihy]; pcy=eico->pcyall[ic][ihy];
		    for(ihx=0; ihx<sf_n(cub->ahx); ihx++) { mcx=eico->mcxall[ic][ihx]; pcx=eico->pcxall[ic][ihx];
			
#ifdef _OPENMP
#pragma omp atomic
#endif
			weop->ihic[ic][iht][ihy][ihx] += 
			    dHIC(weop->swfl[ompith][mcy][mcx],
				 weop->rwfl[ompith][pcy][pcx],
				 wt,
				 meps);
			
		    }
		}
	    }
	    
	} /* cc */

    } else { /* correlation */
	
	for(jc=0;jc<eico->ncofz[iz];jc++) {
	    ic=eico->icofz[iz][jc];
	    
	    for        (iht=0; iht<sf_n(cub->aht); iht++) { wt =eico->tt[iw][iht];
		for    (ihy=0; ihy<sf_n(cub->ahy); ihy++) { mcy=eico->mcyall[ic][ihy]; pcy=eico->pcyall[ic][ihy];
		    for(ihx=0; ihx<sf_n(cub->ahx); ihx++) { mcx=eico->mcxall[ic][ihx]; pcx=eico->pcxall[ic][ihx];
			
#ifdef _OPENMP
#pragma omp atomic
#endif
			weop->ihic[ic][iht][ihy][ihx] += 
			    cHIC(weop->swfl[ompith][mcy][mcx],
				 weop->rwfl[ompith][pcy][pcx],
				 wt);
			
		    }
		}
	    }
	    
	} /* cc */
    }
}

/*------------------------------------------------------------*/
void gensou_inject(sf_complex **wfl,
                   sf_complex **sou,
		   weico3d  eico,
		   int iz,
		   int iw)
/*< inject generalized source into wavefield >*/
{
    int ic, jc, icx, icy;

    for(jc=0;jc<eico->ncofz[iz];jc++) {
        ic=eico->icofz[iz][jc];
        icy = eico->mcyall[ic][0];
        icx = eico->mcxall[ic][0];

#ifdef SF_HAS_COMPLEX_H
        wfl[icy][icx] += sou[iw][ic];
#else
	wfl[icy][icx] = sf_cadd(wfl[icy][icx],sou[iw][ic]);
#endif
    } /* cc */
}

/*------------------------------------------------------------*/
void for_cic(weicub3d cub,
	     sf_complex ***swfl,
	     sf_complex ***rwfl,
	     sf_complex  ***cic,
	     bool conj)
/*< apply forward C.I.C. >*/
{

    int iy, ix, iz;

    for(  iz=sf_n(cub->az)-1; iz>-1; iz--){
	for(  iy=0; iy<sf_n(cub->amy); iy++){
	    for(ix=0; ix<sf_n(cub->amx); ix++){
#ifdef SF_HAS_COMPLEX_H
		if(conj)
		    cic[iz][iy][ix] += cWGH(swfl[iz][iy][ix],rwfl[iz][iy][ix],sf_cmplx(0.0,-0.69));
		else
		    cic[iz][iy][ix] += cAWGH(swfl[iz][iy][ix],rwfl[iz][iy][ix],sf_cmplx(0.0,-0.69));
#else
		if(conj)
		    cic[iz][iy][ix] = sf_cadd(cic[iz][iy][ix],cWGH(swfl[iz][iy][ix],rwfl[iz][iy][ix],sf_cmplx(0.0,-0.69)));
		else
		    cic[iz][iy][ix] = sf_cadd(cic[iz][iy][ix],cAWGH(swfl[iz][iy][ix],rwfl[iz][iy][ix],sf_cmplx(0.0,-0.69)));
#endif
	    } /* loop over ax */
	} /* loop over ay */
    } /* loop over az */
}



/*------------------------------------------------------------*/
void weicic_apply(weiop3d weop,
		  weicub3d cub)
/*< apply CIC >*/
{
    int ix,iy,iz,nz;
    
    nz = sf_n(cub->az);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ix,iy,iz)				\
    shared(weop)
#endif
#ifdef SF_HAS_COMPLEX_H
    CICLOOP(
	;            weop->icic[iz][iy][ix] +=
	crealf(conjf(weop->sone[iz][iy][ix])
	       *     weop->rone[iz][iy][ix]);
	);
#else
    CICLOOP(
	;            weop->icic[iz][iy][ix] +=
	crealf(sf_cmul(conjf(weop->sone[iz][iy][ix]),
		       weop->rone[iz][iy][ix]));
	);
#endif
}

/*------------------------------------------------------------*/
void weieic_apply(int jw,
		  int ibk,
		  weiop3d weop,
		  weico3d eico,
		  weicub3d cub)
/*< apply EIC >*/
{
    int ihx,ihy,ihz,iht,nht,ic,nc;
    sf_complex wt;
    int mcx,pcx;
    int mcy,pcy;
    int mcz,pcz;
    int iw,nw;

    nw=cub->ompnth;  
    iw=ibk*nw+jw;

    nc = sf_n(cub->ac);
    nht = sf_n(cub->aht);

    if(sf_n(cub->ac)>cub->ompnth) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)			\
    private(ic,iht,ihz,ihy,ihx,wt,mcz,pcz,mcy,pcy,mcx,pcx)	\
    shared(cub,eico,weop,iw,nc)
#endif
	for(ic=0; ic<nc; ic++) {
	    if(eico->ccin[ic]) {
		
		for            (iht=0; iht<sf_n(cub->aht); iht++) { wt =eico->tt[iw][iht];
		    for        (ihz=0; ihz<sf_n(cub->ahz); ihz++) { mcz=eico->mczall[ic][ihz]; pcz=eico->pczall[ic][ihz];
			for    (ihy=0; ihy<sf_n(cub->ahy); ihy++) { mcy=eico->mcyall[ic][ihy]; pcy=eico->pcyall[ic][ihy];
			    for(ihx=0; ihx<sf_n(cub->ahx); ihx++) { mcx=eico->mcxall[ic][ihx]; pcx=eico->pcxall[ic][ihx];
				weop->ieic[ic][iht][ihz][ihy][ihx] += 
				    crealf(
					cWGH(weop->sone[mcz][mcy][mcx],weop->rone[pcz][pcy][pcx],wt)
					); 
			    }
			}
		    }
		}
		
	    } /* ccin */
	} /* cc */
	
    } else {
	for(ic=0; ic<sf_n(cub->ac); ic++) {
	    if(eico->ccin[ic]) {
		
#ifdef _OPENMP
#pragma omp parallel for schedule(static)		\
    private(iht,ihz,ihy,ihx,wt,mcz,pcz,mcy,pcy,mcx,pcx) \
    shared(cub,eico,weop,iw,ic)
#endif
		for            (iht=0; iht<nht; iht++) { wt =eico->tt[iw][iht];
		    for        (ihz=0; ihz<sf_n(cub->ahz); ihz++) { mcz=eico->mczall[ic][ihz]; pcz=eico->pczall[ic][ihz];
			for    (ihy=0; ihy<sf_n(cub->ahy); ihy++) { mcy=eico->mcyall[ic][ihy]; pcy=eico->pcyall[ic][ihy];
			    for(ihx=0; ihx<sf_n(cub->ahx); ihx++) { mcx=eico->mcxall[ic][ihx]; pcx=eico->pcxall[ic][ihx];
				weop->ieic[ic][iht][ihz][ihy][ihx] += 
				    crealf(
					cWGH(weop->sone[mcz][mcy][mcx],weop->rone[pcz][pcy][pcx],wt)
					); 
			    }
			}
		    }
		}
		
	    } /* ccin */
	} /* cc */
	
    } /* else */
    
}
