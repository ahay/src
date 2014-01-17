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
#include "weiajs.h"
#include "weiutl.h"
/*^*/

#ifdef SF_HAS_COMPLEX_H
#define cWGH(a,b,c) ((conjf(a)*b)*c)
#define cAWGH(a,b,c) ((a*b)*c)
#else
#define cWGH(a,b,c) sf_cmul(sf_cmul(conjf(a),b),c)
#define cAWGH(a,b,c) sf_cmul(sf_cmul(a,b),c)
#endif

#define  YXLOOP(a) \
    for(iy=0;iy<sf_n(cub->amy);iy++){			\
	for(ix=0;ix<sf_n(cub->amx);ix++){		\
	    {a} }} /*  yx loop */

#define CICLOOP(a) \
    for(iz=0;iz<sf_n(cub->az);iz++){				\
	for(iy=0;iy<sf_n(cub->amy);iy++){			\
	    for(ix=0;ix<sf_n(cub->amx);ix++){			\
		{a} }}} /* zyx loop */

#define HICLOOP(a) \
    for(ic=0;  ic<nc; ic++ ) {					\
	for(iht=0;iht<sf_n(cub->aht);iht++) {			\
	    for(ihy=0;ihy<sf_n(cub->ahy);ihy++) {		\
		for(ihx=0;ihx<sf_n(cub->ahx);ihx++) {		\
		    {a} }}}} /* cc-ht-hy-hx loop */

/*------------------------------------------------------------*/
weiop3d weihicajs_init(weicub3d cub)
/*< initialize HIC mig >*/
{
    weiop3d weop;
    weop = (weiop3d) sf_alloc(1,sizeof(*weop));
    
    weop->wfld = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    weop->ihic = sf_floatalloc4  (sf_n(cub->ahx),sf_n(cub->ahy),sf_n(cub->aht),sf_n(cub->ac));
    weop->ahic = sf_complexalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);

    return weop;
}

/*------------------------------------------------------------*/
void mvahic_inp(weicub3d cub,
                sf_file Feic,
                sf_file Fcoo)
/*< input HIC parameters from files >*/
{
    sf_axis ahx, ahy, ahz, aht, ac;

    ahx = sf_iaxa(Feic,1); sf_setlabel(ahx,"hx");
    ahy = sf_iaxa(Feic,2); sf_setlabel(ahy,"hy");
    ahz = sf_iaxa(Feic,3); sf_setlabel(ahz,"hz");
    aht = sf_iaxa(Feic,4); sf_setlabel(aht,"ht");

    ac = sf_iaxa(Fcoo,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");

    sf_copyaxis(cub->ac,ac);
    sf_copyaxis(cub->ahx,ahx);
    sf_copyaxis(cub->ahy,ahy);
    sf_copyaxis(cub->ahz,ahz);
    sf_copyaxis(cub->aht,aht);
}

/*------------------------------------------------------------*/
void adjsou(weicub3d cub,
            weico3d eico,
            sf_file Feic,
            sf_file Fbwf,
            sf_file Faso,
	    bool conj)
/*< adjoint source construction >*/
{

    int ix, iy, iz, iw, ompith=0;
    sf_complex *****eic, ****bwf, ****asou;

    eic = sf_complexalloc5(sf_n(cub->ahx),sf_n(cub->ahy),sf_n(cub->ahz),sf_n(cub->aht),sf_n(cub->ac));
    bwf = sf_complexalloc4(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az),cub->ompnth);
    asou= sf_complexalloc4(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az),cub->ompnth);

    /*------------------------------------------------------------*/
    /* read extended image */
    sf_complexread(eic[0][0][0][0],sf_n(cub->ahx)*sf_n(cub->ahy)*sf_n(cub->ahz)*sf_n(cub->aht)*sf_n(cub->ac),Feic);

    /*------------------------------------------------------------*/
    /* loop over frequency */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(ompith,iw,iz,iy,ix)		\
    shared(cub,eico)
#endif
    for (iw=0; iw<sf_n(cub->aw); iw++) {
#ifdef _OPENMP
        ompith = omp_get_thread_num();
#pragma omp critical
	if(cub->verb) sf_warning ("eAJS ... (ith=%2d) ... <iw=%3d of %3d>",
				  ompith,iw+1,sf_n(cub->aw));
#endif

	CICLOOP( asou[ompith][iz][iy][ix] = sf_cmplx(0.0,0.0); );

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            sf_seek(Fbwf,sizeof(sf_complex)*cub->nxy*sf_n(cub->az)*iw,SEEK_SET);
            sf_complexread(bwf[ompith][0][0],cub->nxy*sf_n(cub->az),Fbwf);
	}

        adj_eic(cub,eico,bwf[ompith],asou[ompith],eic,conj,iw);

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            /* output wavefield */
            sf_seek(Faso,sizeof(sf_complex)*cub->nxy*sf_n(cub->az)*iw,SEEK_SET);
            sf_complexwrite(asou[ompith][0][0],cub->nxy*sf_n(cub->az),Faso);
        }
    } /* iw */

    free( ****eic); free( ***eic); free( **eic); free( *eic); free( eic);
    free( ***bwf);  free( **bwf);  free( *bwf);  free( bwf);
    free( ***asou); free( **asou); free( *asou); free( asou);
}

/*------------------------------------------------------------*/
void adj_eic(weicub3d cub,
             weico3d eico,
	     sf_complex ***swfl,
	     sf_complex ***asou,
	     sf_complex *****cip,
	     bool cjg,
	     int iw)
/*< apply adjoint E.I.C. >*/
{
    int ihx,ihy,ihz,iht,ic;
    sf_complex wt;
    int mcz,pcz;
    int mcx,pcx;
    int mcy,pcy;
    
    for                (ic= 0; ic <sf_n(cub->ac);  ic++ ) {
        for            (iht=0; iht<sf_n(cub->aht); iht++) { wt =eico->tt[iw][iht];
            for        (ihz=0; ihz<sf_n(cub->ahz); ihz++) { mcz=eico->mczall[ic][ihz]; pcz=eico->pczall[ic][ihz];
                for    (ihy=0; ihy<sf_n(cub->ahy); ihy++) { mcy=eico->mcyall[ic][ihy]; pcy=eico->pcyall[ic][ihy];
                    for(ihx=0; ihx<sf_n(cub->ahx); ihx++) { mcx=eico->mcxall[ic][ihx]; pcx=eico->pcxall[ic][ihx];

#ifdef SF_HAS_COMPLEX_H
			if(cjg)	asou[mcz][mcy][mcx] += cAWGH(swfl[pcz][pcy][pcx],cip[ic][iht][ihz][ihy][ihx],     (wt));
			else    asou[pcz][pcy][pcx] += cAWGH(swfl[mcz][mcy][mcx],cip[ic][iht][ihz][ihy][ihx],conjf(wt));
#else
			if(cjg) asou[mcz][mcy][mcx] = sf_cadd(asou[mcz][mcy][mcx],cAWGH(swfl[pcz][pcy][pcx],cip[ic][iht][ihz][ihy][ihx],     (wt)));
			else    asou[pcz][pcy][pcx] = sf_cadd(asou[pcz][pcy][pcx],cAWGH(swfl[mcz][mcy][mcx],cip[ic][iht][ihz][ihy][ihx],conjf(wt)));
#endif
                    }
                }
            }
        }
    } /* cc */
}

/*------------------------------------------------------------*/
void weihicajs_close(weiop3d weop)
/*< free HIC adjoint source storage >*/
{
    ;                   free(**weop->wfld);free(*weop->wfld);free( weop->wfld);
    free(***weop->ihic);free(**weop->ihic);free(*weop->ihic);free( weop->ihic);
    ;                   free(**weop->ahic);free(*weop->ahic);free( weop->ahic);
}

/*------------------------------------------------------------*/
void weihicajs(weiop3d weop,
	       weicub3d cub,
	       weico3d eico,
	       sf_file Feic,
	       sf_file Fbwf,
	       sf_file Faso,
	       bool conj)
/*< HIC adjoint source construction >*/
{
    int ix,iy,iz,iw,ompith=0;

    /*------------------------------------------------------------*/
    /* read extended image */
    sf_seek(Feic,0,SEEK_SET);
    sf_floatread(weop->ihic[0][0][0],
		 sf_n(cub->ac)*sf_n(cub->aht)*sf_n(cub->ahy)*sf_n(cub->ahx),
		 Feic);

    /*------------------------------------------------------------*/
    /* loop over frequency */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(ompith,iw,iy,ix,iz)		\
    shared(weop,cub,eico)
#endif
    for (iw=0; iw<sf_n(cub->aw); iw++) {
#ifdef _OPENMP
        ompith = omp_get_thread_num();
#pragma omp critical
	if(cub->verb) sf_warning ("hAJS ... (ith=%2d) ... <iw=%3d of %3d>",
				  ompith,iw+1,sf_n(cub->aw));
#endif
	
	/* loop over z */
	for (iz=0; iz<sf_n(cub->az); iz++) {
	    
	    /* zero adjoint source */
	    YXLOOP( weop->ahic[ompith][iy][ix] = sf_cmplx(0.0,0.0); );

#ifdef _OPENMP
#pragma omp critical
#endif
	    {
		/* read wavefield */
		sf_seek(Fbwf,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
		sf_complexread(weop->wfld[ompith][0],cub->nxy,Fbwf);
	    }

	    weihicajs_apply(weop,cub,eico,ompith,iw,iz,conj);
	    
#ifdef _OPENMP
#pragma omp critical
#endif
	    {
		/* output adjoint source */
		sf_seek(Faso,sizeof(sf_complex)*cub->nxy*(sf_n(cub->az)*iw+iz),SEEK_SET);
		sf_complexwrite(weop->ahic[ompith][0],cub->nxy,Faso);
	    }
	    
	} /* iz */
    } /* iw */

}

/*------------------------------------------------------------*/
void weihicajs_apply(weiop3d weop,
		     weicub3d cub,
		     weico3d eico,
		     int ompith,
		     int iw,
		     int iz,
		     bool cjg)
/*< apply adjoint E.I.C. >*/
{
    int ihx,ihy,iht,ic,jc;
    sf_complex wt;
    int mcx,pcx;
    int mcy,pcy;

    if(cjg) {

	for(jc=0;jc<eico->ncofz[iz];jc++) {
	    ic=eico->icofz[iz][jc];

	    for        (iht=0; iht<sf_n(cub->aht); iht++) { wt =eico->tt[iw][iht];
		for    (ihy=0; ihy<sf_n(cub->ahy); ihy++) { mcy=eico->mcyall[ic][ihy]; pcy=eico->pcyall[ic][ihy];
		    for(ihx=0; ihx<sf_n(cub->ahx); ihx++) { mcx=eico->mcxall[ic][ihx]; pcx=eico->pcxall[ic][ihx];

#ifdef SF_HAS_COMPLEX_H
			weop->ahic[ompith][mcy][mcx] += cAWGH(weop->wfld[ompith][pcy][pcx],
							      weop->ihic[ic][iht][ihy][ihx],     (wt));
#else
			weop->ahic[ompith][mcy][mcx] = sf_cadd(weop->ahic[ompith][mcy][mcx],
							       cAWGH(weop->wfld[ompith][pcy][pcx],
								     weop->ihic[ic][iht][ihy][ihx],     (wt)));
#endif
		    }
		}
	    }
	}
	
    } else {
 	
	for(jc=0;jc<eico->ncofz[iz];jc++) {
	    ic=eico->icofz[iz][jc];
	    
	    for        (iht=0; iht<sf_n(cub->aht); iht++) { wt =eico->tt[iw][iht];
		for    (ihy=0; ihy<sf_n(cub->ahy); ihy++) { mcy=eico->mcyall[ic][ihy]; pcy=eico->pcyall[ic][ihy];
		    for(ihx=0; ihx<sf_n(cub->ahx); ihx++) { mcx=eico->mcxall[ic][ihx]; pcx=eico->pcxall[ic][ihx];
#ifdef SF_HAS_COMPLEX_H
			weop->ahic[ompith][pcy][pcx] += cAWGH(weop->wfld[ompith][mcy][mcx],
							      weop->ihic[ic][iht][ihy][ihx],conjf(wt));
#else
			weop->ahic[ompith][pcy][pcx] = sf_cadd(weop->ahic[ompith][pcy][pcx],
							       cAWGH(weop->wfld[ompith][mcy][mcx],
								     weop->ihic[ic][iht][ihy][ihx],conjf(wt)));
#endif
		    }
		}
	    }
	}

    } /* end else */   
}

