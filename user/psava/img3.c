/* Imaging condition for shot-profile migration */

/*
  Copyright (C) 2007 Colorado School of Mines
  
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

#include "img3.h"

#include "weutil.h"
/*^*/

#define  LOOP(a) for( imy=0; imy<cub->amy.n; imy++){ \
                 for( imx=0; imx<cub->amx.n; imx++){ \
		     {a} \
		 }}

#define XLOOP(a) for( imz=0; imz<cub->amz.n; imz++){ \
                 for( imy=0; imy<cub->amy.n; imy++){ \
                 for( imx=0; imx<cub->amx.n; imx++){ \
		     {a} \
		 }}}

#define HLOOP(a) for( ihz=img->LOz; ihz<img->HIz; ihz++){ \
                 for( ihy=img->LOy; ihy<img->HIy; ihy++){ \
                 for( ihx=img->LOx; ihx<img->HIx; ihx++){ \
		     {a} \
		 }}}

#define CLOOP(a) for( icz=0; icz<img->acz.n; icz++){ \
                 for( icy=0; icy<img->acy.n; icy++){ \
                 for( icx=0; icx<img->acx.n; icx++){ \
		     {a} \
		 }}}

#define IND(ihx,ihy,ihz) \
    (ihz-img->LOz)*(img->ahx.n*img->ahy.n) +  \
    (ihy-img->LOy)* img->ahx.n             +  \
    (ihx-img->LOx)

#define EICIND(ihx,ihy,ihz,iht) \
    iht * (img->ahx.n*img->ahy.n*img->ahz.n) + \
    (ihz-img->LOz)*(img->ahx.n*img->ahy.n)   + \
    (ihy-img->LOy)* img->ahx.n               + \
    (ihx-img->LOx)


#define MM(i,a) (SF_MIN(SF_MAX(i,0),a.n-1))

#define INBOUND(imin,imax,i) ((i>=imin && i<imax)?true:false)

/*------------------------------------------------------------*/
static float corr(sf_complex a, sf_complex b)
{
    sf_complex c;
#ifdef SF_HAS_COMPLEX_H
    c = conjf(a)*b;
#else
    c = sf_cmul(conjf(a),b);
#endif
    return crealf(c);
}

/*------------------------------------------------------------*/
static float wcorr(sf_complex a, sf_complex b, sf_complex w)
{
    sf_complex c;
#ifdef SF_HAS_COMPLEX_H
    c = conjf(a)*b*w;
#else
    c = sf_cmul(sf_cmul(conjf(a),b),w);
#endif
    return crealf(c);
}

/*------------------------------------------------------------*/
img3d img3_init(cub3d cub,
		sf_fslice imag,
		sf_fslice cigs,
		int jcx_,
		int jcy_,
		int jcz_
    )
/*< initialize I.C. >*/
{
    int imx,imy,imz;
    int icx,icy,icz;
    int ompith;
    
    /*------------------------------------------------------------*/
    img3d img;
    img = (img3d) sf_alloc(1,sizeof(*img));
    
    img->cigs = cigs;
    img->imag = imag;
    
    img->jcx = jcx_;
    img->jcy = jcy_;
    img->jcz = jcz_;
    
    img->acx.n = cub->amx.n / img->jcx;
    img->acy.n = cub->amy.n / img->jcy;
    img->acz.n = cub->amz.n / img->jcz;
    
    /* allocate wavefield storage */
    img->qs = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->amz.n,cub->ompnth);
    img->qr = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->amz.n,cub->ompnth);
    
    /* allocate image storage */
    img->qi = sf_floatalloc4(cub->amx.n,cub->amy.n,cub->amz.n,cub->ompnth);
    
    for(ompith=0; ompith<cub->ompnth; ompith++){
	XLOOP( img->qs[ompith][imz][imy][imx] = sf_cmplx(0.0,0.0); 
	       img->qr[ompith][imz][imy][imx] = sf_cmplx(0.0,0.0); 
	       img->qi[ompith][imz][imy][imx] = 0.0;
	    );
    }

    /* allocate cigs storage */
    img->qc = sf_floatalloc4(img->acx.n,img->acy.n,img->acz.n,cub->ompnth);
    img->qt = sf_floatalloc4(img->acx.n,img->acy.n,img->acz.n,cub->ompnth);
    for(ompith=0; ompith<cub->ompnth; ompith++){
	CLOOP( img->qc[ompith][icz][icy][icx] = 0.0; );
	CLOOP( img->qt[ompith][icz][icy][icx] = 0.0; );
    }

    return img;
}

/*------------------------------------------------------------*/
int cipmin(int nn,
	   int jc,
	   int nc,
	   int nh)
/*< find min eCIP index >*/
{
    int ic;
    for(ic=nc;ic*jc-nh>0;ic--) {
	;
    }
    ic++;
/*    i = (1.0*nh/jc)+1;*/
    return ic;
}

/*------------------------------------------------------------*/
int cipmax(int nn,
	   int jc,
	   int nc,
	   int nh)
/*< find max eCIP index >*/
{
    int ic;
    for(ic=0;ic*jc+nh<nn-1;ic++) {
	;
    }
    ic--;
/*    i = (1.0*(nn-1-nh)/jc)+1;*/
    return ic;
}

/*------------------------------------------------------------*/
img3d img3o_init(cub3d cub,
		 sf_fslice imag,
		 sf_fslice cigs,
		 int jcx_,
		 int jcy_,
		 int jcz_
    )
/*< initialize zero-lag I.C. >*/
{
    img3d img;
    img=img3_init(cub,imag,cigs,jcx_,jcy_,jcz_);

    sf_fslice_put(img->cigs,0,img->qc[0][0][0]);

    return img;
}

/*------------------------------------------------------------*/
img3d img3x_init(cub3d cub,
		 sf_fslice imag,
		 sf_fslice cigs,
		 int jcx_,
		 int jcy_,
		 int jcz_,
		 sf_axis ahx_,
		 sf_axis ahy_,
		 sf_axis ahz_
    )
/*< initialize x-lag I.C. >*/
{
    int ihx,ihy,ihz;
    img3d img;
    int nhmax;
    img=img3_init(cub,imag,cigs,jcx_,jcy_,jcz_);

    img->ahx = sf_nod(ahx_);
    img->ahy = sf_nod(ahy_);
    img->ahz = sf_nod(ahz_);

    img->LOx = floor(img->ahx.o/img->ahx.d); img->HIx = img->LOx + img->ahx.n;
    img->LOy = floor(img->ahy.o/img->ahy.d); img->HIy = img->LOy + img->ahy.n;
    img->LOz = floor(img->ahz.o/img->ahz.d); img->HIz = img->LOz + img->ahz.n;

    nhmax = SF_MAX(SF_ABS(img->LOx),SF_ABS(img->HIx));
    img->CLOx = cipmin(cub->amx.n,img->jcx,img->acx.n,nhmax); img->CLOx = MM(img->CLOx,img->acx);
    img->CHIx = cipmax(cub->amx.n,img->jcx,img->acx.n,nhmax); img->CHIx = MM(img->CHIx,img->acx);

    nhmax = SF_MAX(SF_ABS(img->LOy),SF_ABS(img->HIy));    
    img->CLOy = cipmin(cub->amy.n,img->jcy,img->acy.n,nhmax); img->CLOy = MM(img->CLOy,img->acy);
    img->CHIy = cipmax(cub->amy.n,img->jcy,img->acy.n,nhmax); img->CHIy = MM(img->CHIy,img->acy);

    nhmax = SF_MAX(SF_ABS(img->LOz),SF_ABS(img->HIz));
    img->CLOz = cipmin(cub->amz.n,img->jcz,img->acz.n,nhmax); img->CLOz = MM(img->CLOz,img->acz);
    img->CHIz = cipmax(cub->amz.n,img->jcz,img->acz.n,nhmax); img->CHIz = MM(img->CHIz,img->acz);

    img->CHIx = SF_MAX(1,img->CHIx);
    img->CHIy = SF_MAX(1,img->CHIy);
    img->CHIz = SF_MAX(1,img->CHIz);

    if(cub->verb) {
	sf_warning("CLOx=%3d CHIx=%3d",img->CLOx,img->CHIx);
	sf_warning("CLOy=%3d CHIy=%3d",img->CLOy,img->CHIy);
	sf_warning("CLOz=%3d CHIz=%3d",img->CLOz,img->CHIz);
    }

    HLOOP(
	sf_fslice_put(img->cigs,IND(ihx,ihy,ihz),img->qc[0][0][0]);
	);

    return img;
}

/*------------------------------------------------------------*/
img3d img3t_init(cub3d cub,
		 sf_fslice imag,
		 sf_fslice cigs,
		 int jcx_,
		 int jcy_,
		 int jcz_,
		 sf_axis aht_
    )
/*< initialize t-lag I.C. >*/
{
    int  iht,iw;
    float ht, w;
    img3d img;
    img=img3_init(cub,imag,cigs,jcx_,jcy_,jcz_);

    img->aht = sf_nod(aht_);

    for( iht=0; iht<img->aht.n; iht++) {
	sf_fslice_put(img->cigs,iht,img->qc[0][0][0]);
    }

    /* precompute phase shift */
    img->tt = sf_complexalloc2(img->aht.n,cub->aw.n);
    for (iw=0; iw<cub->aw.n; iw++) {
	w = cub->aw.o+iw*cub->aw.d;
	for (iht=0; iht<img->aht.n; iht++) {
	    ht = img->aht.o+iht*img->aht.d;
	    ht *= -2*w;
	    img->tt[iw][iht] = sf_cmplx(cosf(ht),sinf(ht));
	}
    }

    return img;
}

/*------------------------------------------------------------*/
img3d img3e_init(cub3d cub,
		 sf_fslice imag,
		 sf_fslice cigs,
		 int jcx_,
		 int jcy_,
		 int jcz_,
		 sf_axis ahx_,
		 sf_axis ahy_,
		 sf_axis ahz_,
		 sf_axis aht_
    )
/*< initialize E.I.C. >*/
{
    int ihx,ihy,ihz,iht,iw;
    float ht, w;
    img3d img;
    int nhmax;
    img=img3_init(cub,imag,cigs,jcx_,jcy_,jcz_);

    img->ahx = sf_nod(ahx_);
    img->ahy = sf_nod(ahy_);
    img->ahz = sf_nod(ahz_);
    img->aht = sf_nod(aht_);

    img->LOx = floor(img->ahx.o/img->ahx.d); img->HIx = img->LOx + img->ahx.n;
    img->LOy = floor(img->ahy.o/img->ahy.d); img->HIy = img->LOy + img->ahy.n;
    img->LOz = floor(img->ahz.o/img->ahz.d); img->HIz = img->LOz + img->ahz.n;
    
    nhmax = SF_MAX(SF_ABS(img->LOx),SF_ABS(img->HIx));
    img->CLOx = cipmin(cub->amx.n,img->jcx,img->acx.n,nhmax); img->CLOx = MM(img->CLOx,img->acx);
    img->CHIx = cipmax(cub->amx.n,img->jcx,img->acx.n,nhmax); img->CHIx = MM(img->CHIx,img->acx);

    nhmax = SF_MAX(SF_ABS(img->LOy),SF_ABS(img->HIy));    
    img->CLOy = cipmin(cub->amy.n,img->jcy,img->acy.n,nhmax); img->CLOy = MM(img->CLOy,img->acy);
    img->CHIy = cipmax(cub->amy.n,img->jcy,img->acy.n,nhmax); img->CHIy = MM(img->CHIy,img->acy);

    nhmax = SF_MAX(SF_ABS(img->LOz),SF_ABS(img->HIz));
    img->CLOz = cipmin(cub->amz.n,img->jcz,img->acz.n,nhmax); img->CLOz = MM(img->CLOz,img->acz);
    img->CHIz = cipmax(cub->amz.n,img->jcz,img->acz.n,nhmax); img->CHIz = MM(img->CHIz,img->acz);

    img->CHIx = SF_MAX(1,img->CHIx);
    img->CHIy = SF_MAX(1,img->CHIy);
    img->CHIz = SF_MAX(1,img->CHIz);

    if(cub->verb) {
	sf_warning("CLOx=%3d CHIx=%3d",img->CLOx,img->CHIx);
	sf_warning("CLOy=%3d CHIy=%3d",img->CLOy,img->CHIy);
	sf_warning("CLOz=%3d CHIz=%3d",img->CLOz,img->CHIz);
    }

    for( iht=0; iht<img->aht.n; iht++) {
	HLOOP(
	    sf_fslice_put(img->cigs,EICIND(ihx,ihy,ihz,iht),img->qc[0][0][0]);
	    );
    }
    
    /* precompute phase shift */
    img->tt = sf_complexalloc2(img->aht.n,cub->aw.n);
    for (iw=0; iw<cub->aw.n; iw++) {
	w = cub->aw.o+iw*cub->aw.d;
	for (iht=0; iht<img->aht.n; iht++) {
	    ht = img->aht.o+iht*img->aht.d;
	    ht *= -2*w;
	    img->tt[iw][iht] = sf_cmplx(cosf(ht),sinf(ht));
	}
    }

    return img;
}

/*------------------------------------------------------------*/
img3d img3h_init(cub3d cub,
		 sf_fslice imag,
		 sf_fslice cigs,
		 int jcx_,
		 int jcy_,
		 int jcz_,
		 sf_axis ahh_,
		 sf_axis aha_,
		 sf_axis ahb_,
		 float vpvs_
    )
/*< initialize abs x-lag I.C. >*/
{
    int  ihh;
    img3d img;
    img=img3_init(cub,imag,cigs,jcx_,jcy_,jcz_);

    img->ahh = sf_nod(ahh_);
    img->aha = sf_nod(aha_);
    img->ahb = sf_nod(ahb_);

    img->vpvs = vpvs_;

    for( ihh=0; ihh<img->ahh.n; ihh++) {
	sf_fslice_put(img->cigs,ihh,img->qc[0][0][0]);
    }

    return img;
}

/*------------------------------------------------------------*/
void img3_cout(img3d img,
	       int i,
	       int ompith)
/*< update cigs on disk >*/
{
    int icx,icy,icz;

#ifdef _OPENMP	    
#pragma omp critical
#endif
    {
	sf_fslice_get(img->cigs,i,img->qt[ompith][0][0]);
	CLOOP(
	    img->qt[ompith][icz][icy][icx] +=
	    img->qc[ompith][icz][icy][icx];
	    );
	sf_fslice_put(img->cigs,i,img->qt[ompith][0][0]);
    }   
}

/*------------------------------------------------------------*/
void img3o(cub3d cub,
	   img3d img,
	   int iw,
	   int ompith)
/*< apply zero-lag I.C. >*/
{
    int imx,imy,imz;
    int icx,icy,icz;

    /* imag */
    XLOOP(;    img->qi[ompith][imz][imy][imx] +=
	  corr(img->qs[ompith][imz][imy][imx],    
	       img->qr[ompith][imz][imy][imx]);
	);

    CLOOP(;    img->qc[ompith][         icz][         icy][         icx] =
	  corr(img->qs[ompith][icz*img->jcz][icy*img->jcy][icx*img->jcx], 
	       img->qr[ompith][icz*img->jcz][icy*img->jcy][icx*img->jcx]);
	);

    /* cigs */
    img3_cout(img,0,ompith);
}

/*------------------------------------------------------------*/
void oldimg3x(cub3d cub,
	   img3d img,
	   int iw,
	   int ompith)
/*< apply x-lag I.C. >*/
{
    int imx, imy, imz;
    int icx, icy, icz;
    int ihx, ihy, ihz;
    int imys,imyr,imzs;
    int imxs,imxr,imzr;

    /* imag */
    XLOOP(;    img->qi[ompith][imz][imy][imx] +=
	  corr(img->qs[ompith][imz][imy][imx], 
	       img->qr[ompith][imz][imy][imx]);
	);

    for        ( ihz=img->LOz; ihz<img->HIz; ihz++){ 
	for    ( ihy=img->LOy; ihy<img->HIy; ihy++){
	    for( ihx=img->LOx; ihx<img->HIx; ihx++){
		
		for( icz=0; icz<img->acz.n; icz++){ 
		    imzs = icz*img->jcz - ihz;
		    imzr = icz*img->jcz + ihz;
		    if ( INBOUND(0,cub->amz.n,imzs) &&
			 INBOUND(0,cub->amz.n,imzr)) {		    

			for( icy=0; icy<img->acy.n; icy++){
			    imys = icy*img->jcy - ihy;
			    imyr = icy*img->jcy + ihy;
			    if ( INBOUND(0,cub->amy.n,imys) &&
				 INBOUND(0,cub->amy.n,imyr)) {
				
				for( icx=0; icx<img->acx.n; icx++){
				    imxs = icx*img->jcx - ihx;
				    imxr = icx*img->jcx + ihx;
				    if ( INBOUND(0,cub->amx.n,imxs) &&
					 INBOUND(0,cub->amx.n,imxr)) {
					
					;        img->qc[ompith][icz] [icy] [icx] =
					    corr(img->qs[ompith][imzs][imys][imxs], 
						 img->qr[ompith][imzr][imyr][imxr]);
				    }
				} /* cx */
			    }
			}         /* cy */
		    }
		}                 /* cz */
		
		img3_cout(img,IND(ihx,ihy,ihz),ompith);

	    } /* hx */
	}     /* hy */
    }         /* hz */
}

/*------------------------------------------------------------*/
void img3x(cub3d cub,
	   img3d img,
	   int iw,
	   int ompith)
/*< apply x-lag I.C. >*/
{
    int imx, imy, imz;
    int icx, icy, icz;
    int ihx, ihy, ihz;
    int imys,imyr,imzs;
    int imxs,imxr,imzr;

    /* imag */
    XLOOP(;    img->qi[ompith][imz][imy][imx] +=
	  corr(img->qs[ompith][imz][imy][imx], 
	       img->qr[ompith][imz][imy][imx]);
	);

    for        ( ihz=img->LOz; ihz<img->HIz; ihz++){ 
	for    ( ihy=img->LOy; ihy<img->HIy; ihy++){
	    for( ihx=img->LOx; ihx<img->HIx; ihx++){
		
		for( icz=img->CLOz; icz<=img->CHIz; icz++){ 
		    imzs = icz*img->jcz - ihz;
		    imzr = icz*img->jcz + ihz;

			for( icy=img->CLOy; icy<img->CHIy; icy++){
			    imys = icy*img->jcy - ihy;
			    imyr = icy*img->jcy + ihy;
				
				for( icx=img->CLOx; icx<=img->CHIx; icx++){
				    imxs = icx*img->jcx - ihx;
				    imxr = icx*img->jcx + ihx;

					;        img->qc[ompith][icz] [icy] [icx] =
					    corr(img->qs[ompith][imzs][imys][imxs], 
						 img->qr[ompith][imzr][imyr][imxr]);
				} /* cx */
			}         /* cy */
		}                 /* cz */
		
		img3_cout(img,IND(ihx,ihy,ihz),ompith);

	    } /* hx */
	}     /* hy */
    }         /* hz */
}

/*------------------------------------------------------------*/
void img3t(cub3d cub,
	   img3d img,
	   int iw,
	   int ompith)
/*< apply t-lag I.C. >*/
{
    int imx,imy,imz,iht;
    int icx,icy,icz;
    sf_complex wt;

    /* imag */
    XLOOP(;    img->qi[ompith][imz][imy][imx] +=
	  corr(img->qs[ompith][imz][imy][imx], 
	       img->qr[ompith][imz][imy][imx]);
	);

    /* cigs */
    for(iht=0; iht<img->aht.n; iht++) {
	wt = img->tt[iw][iht];
	
	CLOOP(;     img->qc[ompith][icz         ][icy         ][icx         ] = 
	      wcorr(img->qs[ompith][icz*img->jcz][icy*img->jcy][icx*img->jcx], 
		    img->qr[ompith][icz*img->jcz][icy*img->jcy][icx*img->jcx],
		    wt);
	    ); 
	img3_cout(img,iht,ompith);
	
    } /* iht */
}

/*------------------------------------------------------------*/
void oldimg3e(cub3d cub,
	   img3d img,
	   int iw,
	   int ompith)
/*< apply E.I.C. >*/
{
    int imx,imy,imz,iht;
    int icx,icy,icz;

    int ihx, ihy, ihz;
    int imys,imyr,imzs;
    int imxs,imxr,imzr;

    sf_complex wt;

    /* imag */
    XLOOP(;    img->qi[ompith][imz][imy][imx] +=
	  corr(img->qs[ompith][imz][imy][imx], 
	       img->qr[ompith][imz][imy][imx]);
	);
    
    for (iht=0; iht<img->aht.n; iht++) {
	wt = img->tt[iw][iht];

	for        ( ihz=img->LOz; ihz<img->HIz; ihz++){ 
	    for    ( ihy=img->LOy; ihy<img->HIy; ihy++){
		for( ihx=img->LOx; ihx<img->HIx; ihx++){
		    
		    for( icz=0; icz<img->acz.n; icz++){ 
			imzs = icz*img->jcz - ihz;
			imzr = icz*img->jcz + ihz;
			if ( INBOUND(0,cub->amz.n,imzs) &&
			     INBOUND(0,cub->amz.n,imzr)) {		    
			    
			    for( icy=0; icy<img->acy.n; icy++){
				imys = icy*img->jcy - ihy;
				imyr = icy*img->jcy + ihy;
				if ( INBOUND(0,cub->amy.n,imys) &&
				     INBOUND(0,cub->amy.n,imyr)) {
				    
				    for( icx=0; icx<img->acx.n; icx++){
					imxs = icx*img->jcx - ihx;
					imxr = icx*img->jcx + ihx;
					if ( INBOUND(0,cub->amx.n,imxs) &&
					     INBOUND(0,cub->amx.n,imxr)) {
					    
					    img->qc[ompith][icz] [icy] [icx] =
						wcorr(img->qs[ompith][imzs][imys][imxs], 
						      img->qr[ompith][imzr][imyr][imxr],
						      wt);
					}
				    } /* cx */
				}
			    }         /* cy */
			}
		    }                 /* cz */
		    
		    img3_cout(img,EICIND(ihx,ihy,ihz,iht),ompith);
		    
		} /* hx */
	    }     /* hy */
	}         /* hz */
    }             /* ht */
}

/*------------------------------------------------------------*/
void img3e(cub3d cub,
	   img3d img,
	   int iw,
	   int ompith)
/*< apply E.I.C. >*/
{
    int imx,imy,imz,iht;
    int icx,icy,icz;

    int ihx, ihy, ihz;
    int imys,imyr,imzs;
    int imxs,imxr,imzr;

    sf_complex wt;

    /* imag */
    XLOOP(;    img->qi[ompith][imz][imy][imx] +=
	  corr(img->qs[ompith][imz][imy][imx], 
	       img->qr[ompith][imz][imy][imx]);
	);
    
    for (iht=0; iht<img->aht.n; iht++) {
	wt = img->tt[iw][iht];

	for        ( ihz=img->LOz; ihz<img->HIz; ihz++){ 
	    for    ( ihy=img->LOy; ihy<img->HIy; ihy++){
		for( ihx=img->LOx; ihx<img->HIx; ihx++){
		    
		    for( icz=img->CLOz; icz<=img->CHIz; icz++){ 
			imzs = icz*img->jcz - ihz;
			imzr = icz*img->jcz + ihz;
			    
			    for( icy=img->CLOy; icy<img->CHIy; icy++){
				imys = icy*img->jcy - ihy;
				imyr = icy*img->jcy + ihy;
				    
				for( icx=img->CLOx; icx<=img->CHIx; icx++){
					imxs = icx*img->jcx - ihx;
					imxr = icx*img->jcx + ihx;
					    
					    img->qc[ompith][icz] [icy] [icx] =
						wcorr(img->qs[ompith][imzs][imys][imxs], 
						      img->qr[ompith][imzr][imyr][imxr],
						      wt);
				    } /* cx */
			    }         /* cy */
		    }                 /* cz */
		    
		    img3_cout(img,EICIND(ihx,ihy,ihz,iht),ompith);
		    
		} /* hx */
	    }     /* hy */
	}         /* hz */
    }             /* ht */
}

/*------------------------------------------------------------*/
void img3h(cub3d cub,
	   img3d img,
	   int iw,
	   int ompith)
/*< apply abs-lag I.C. >*/
{
    int imx,imy,imz;
    int ihh,iha,ihb;
    int icx,icy,icz;

    int dsx,dsy,dsz;
    int drx,dry,drz;

    int isx,isy,isz;
    int irx,iry,irz;

    float hh,aa,bb; /* aa,bb in radians */
    float hx,hy,hz;
    sf_complex cs,cr;
    float hscale;

    /* imag */
    XLOOP(;    img->qi[ompith][imz][imy][imx] +=
	  corr(img->qs[ompith][imz][imy][imx],
	       img->qr[ompith][imz][imy][imx] );
	);
    
    /* cigs */
    for(ihh=0; ihh<img->ahh.n; ihh++) {                /* absolute lag */
	hh = img->ahh.o + ihh * img->ahh.d;
	
	CLOOP( img->qc[ompith][icz][icy][icx] =0; );

	hscale=1.;
	if(ihh>0) hscale=(hh+img->ahh.d)*img->aha.d;
	
	for(ihb=0; ihb<img->ahb.n; ihb++) {        /* latitude  */
	    bb = img->ahb.o + ihb * img->ahb.d;
	    for(iha=0; iha<img->aha.n; iha++) {    /* longitude */
		aa = img->aha.o + iha * img->aha.d;
		
		hz = hh * sin(aa);
		hx = hh * cos(aa) * cos(bb);
		hy = hh * cos(aa) * sin(bb);
		
		/* nearest neighbour - source dh */
		dsx = (int)( (2.*img->vpvs/(1.+img->vpvs))*hx / cub->amx.d);
		dsy = (int)( (2.*img->vpvs/(1.+img->vpvs))*hy / cub->amy.d);
		dsz = (int)( (2.*img->vpvs/(1.+img->vpvs))*hz / cub->amz.d);
		
		/* nearest neighbour - receiver dh */
		drx = (int)( (2./(1.+img->vpvs))*hx / cub->amx.d);
		dry = (int)( (2./(1.+img->vpvs))*hy / cub->amy.d);
		drz = (int)( (2./(1.+img->vpvs))*hz / cub->amz.d);
		
		CLOOP(
		    isx=MM(icx*img->jcx-dsx,cub->amx); 
		    isy=MM(icy*img->jcy-dsy,cub->amy);
		    isz=MM(icz*img->jcz-dsz,cub->amz);
		    cs = img->qs[ompith][isz][isy][isx];
		    
		    irx=MM(icx*img->jcx+drx,cub->amx);
		    iry=MM(icy*img->jcy+dry,cub->amy);
		    irz=MM(icz*img->jcz+drz,cub->amz);
		    cr = img->qr[ompith][irz][iry][irx];
		    
		    /* 
		       cs = qs @ x-hx,y-hy,z-hz
		       cr = qr @ x+hx,y+hy,z+hz
		    */
		    img->qc[ompith][icz][icy][icx] += hscale * corr(cs,cr);  
		    ); /* cx,cy,cz */
	    } /* aa */
	}     /* bb */	
	img3_cout(img,ihh,ompith);

    }         /* hh */
}

/*------------------------------------------------------------*/
void img3o_close(cub3d cub,
		 img3d img,
		 sf_fslice imag,
		 sf_fslice cigs)
/*< deallocate zero-lag I.C. >*/
{
    img3_close(cub,img);
}

/*------------------------------------------------------------*/
void img3t_close(cub3d cub,
		 img3d img,
		 sf_fslice imag,
		 sf_fslice cigs)
/*< deallocate t-lag I.C. >*/
{
    img3_close(cub,img);

    free(*img->tt); 
    free( img->tt);
}

/*------------------------------------------------------------*/
void img3h_close(cub3d cub,
		 img3d img,
		 sf_fslice imag,
		 sf_fslice cigs)
/*< deallocate abs-lag I.C. >*/
{
    img3_close(cub,img);
}

/*------------------------------------------------------------*/
void img3x_close(cub3d cub,
		 img3d img,
		 sf_fslice imag,
		 sf_fslice cigs)
/*< deallocate x-lag I.C. >*/
{
    img3_close(cub,img);
}

/*------------------------------------------------------------*/
void img3e_close(cub3d cub,
		 img3d img,
		 sf_fslice imag,
		 sf_fslice cigs)
/*< deallocate x-lag I.C. >*/
{
    img3_close(cub,img);

    free(*img->tt); 
    free( img->tt);
}

/*------------------------------------------------------------*/
void img3_close(cub3d cub,
		img3d img)
/*< deallocate >*/
{
    int imx,imy,imz;
    int ompith;
    
    if(cub->ompnth>0) {
	for(ompith=1; ompith<cub->ompnth; ompith++) {
	    XLOOP( img->qi[0][imz][imy][imx] +=
		   img->qi[ompith][imz][imy][imx]; );
	}
    }
    sf_fslice_put(img->imag,0,img->qi[0][0][0]);

    free(***img->qt); free(**img->qt); free(*img->qt); free(img->qt);
    free(***img->qc); free(**img->qc); free(*img->qc); free(img->qc);
    free(***img->qs); free(**img->qs); free(*img->qs); free(img->qs);
    free(***img->qr); free(**img->qr); free(*img->qr); free(img->qr);
    free(***img->qi) ;free(**img->qi); free(*img->qi); free(img->qi);
}

/*------------------------------------------------------------*/
void img3store(cub3d cub,
	       img3d img,
	       int imz,
	       sf_complex ***ww_s,
	       sf_complex ***ww_r,
	       int ompith
    )
/*< store wavefield >*/
{
    int imx,imy;

    LOOP( img->qs[ompith][imz][imy][imx] = ww_s[ompith][imy][imx];
	  img->qr[ompith][imz][imy][imx] = ww_r[ompith][imy][imx]; );
}

