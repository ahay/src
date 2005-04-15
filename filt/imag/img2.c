/* Imaging condition for shot-profile migration */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "img2.h"

#include "slice.h"
/*^*/

#define LOOP(a) for( imy=0; imy< amy.n; imy++){ \
                for( imx=0; imx< amx.n; imx++){ {a} }}

#define MLOOP(a) for( imz=0; imz< amz.n; imz++){ \
                 for( imy=0; imy< amy.n; imy++){ \
                 for( imx=0; imx< amx.n; imx++){ {a} }}}

#define HLOOP(a) for( ihz=LOz; ihz<HIz; ihz++){ \
                 for( ihy=LOy; ihy<HIy; ihy++){ \
                 for( ihx=LOx; ihx<HIx; ihx++){ {a} }}}

#define CLOOP(a) for( icz=0; icz< acz.n; icz++){ \
                 for( icy=0; icy< acy.n; icy++){ \
                 for( icx=0; icx< acx.n; icx++){ {a} }}}

#define XLOOP(a) for(imz  = abs(ihz); imz<amz.n-abs(ihz) ; imz++){ \
                     imzs = imz - ihz; \
                     imzr = imz + ihz; \
                 for(imy  = abs(ihy); imy<amy.n-abs(ihy) ; imy++){ \
                     imys = imy - ihy; \
                     imyr = imy + ihy; \
                 for(imx  = abs(ihx); imx<amx.n-abs(ihx) ; imx++){ \
                     imxs = imx - ihx; \
                     imxr = imx + ihx; \
                    {a} \
                 }}}

#define IND(ihx,ihy,ihz) (ihz-LOz)*(ahx.n*ahy.n) + \
                         (ihy-LOy)* ahx.n        + \
                         (ihx-LOx);

#define MM(i,a) SF_MIN(SF_MAX(i,0),a.n-1)

static axa amx,amy,amz;
static axa acx,acy,acz;
static int jcx,jcy,jcz;
static axa ahx,ahy,ahz;
static axa aht;
static axa ahh,aha,ahb;
static axa aw;

static float complex  **tt; /* phase shift for time offset */
static float complex ***qs,***qr;
static float         ***qi;

static float         ***qo;
static float      ******qx;
static float        ****qt;
static float        ****qh;

static int LOx,HIx;
static int LOy,HIy;
static int LOz,HIz;

/*------------------------------------------------------------*/

void img2store( int imz,
		float complex **ww_s,
		float complex **ww_r
    )
/*< store wavefield >*/
{
    int imx,imy;

    LOOP( qs[imz][imy][imx] = ww_s[imy][imx];
	  qr[imz][imy][imx] = ww_r[imy][imx]; );
}

/*------------------------------------------------------------*/
void img2o_init(axa amx_,
		axa amy_,
		axa amz_,
		int jcx_,
		int jcy_,
		int jcz_,
		fslice imag
    )
/*< initialize o-offset imaging condition >*/
{
    int imx,imy,imz;
    int icx,icy,icz;

    amx = amx_;
    amy = amy_;
    amz = amz_;

    jcx = jcx_; acx.n = amx.n / jcx;
    jcy = jcy_; acy.n = amy.n / jcy;
    jcz = jcz_; acz.n = amz.n / jcz;

    /* allocate wavefield storage */
    qs = sf_complexalloc3(amx.n,amy.n,amz.n);
    qr = sf_complexalloc3(amx.n,amy.n,amz.n);
    MLOOP( qs[imz][imy][imx] = 0.0; 
	   qr[imz][imy][imx] = 0.0; );

    /* allocate image storage */
    qi = sf_floatalloc3(amx.n,amy.n,amz.n);
    MLOOP( qi[imz][imy][imx] = 0.0; );

    /* allocate cigs storage */
    qo = sf_floatalloc3(acx.n,acy.n,acz.n);
    CLOOP( qo[icz][icy][icx] = 0.0; );
}

void img2x_init(axa amx_,
		axa amy_,
		axa amz_,
		int jcx_,
		int jcy_,
		int jcz_,
		axa ahx_,
		axa ahy_,
		axa ahz_,
		fslice imag
    )
/*< initialize x-offset imaging condition >*/
{
    int imx,imy,imz;
    int icx,icy,icz;
    int ihx,ihy,ihz;

    amx = amx_;
    amy = amy_;
    amz = amz_;

    jcx = jcx_; acx.n = amx.n / jcx;
    jcy = jcy_; acy.n = amy.n / jcy;
    jcz = jcz_; acz.n = amz.n / jcz;
    
    ahx = ahx_;
    ahy = ahy_;
    ahz = ahz_;

    LOx = floor(ahx.o/ahx.d); HIx = LOx + ahx.n;
    LOy = floor(ahy.o/ahy.d); HIy = LOy + ahy.n;
    LOz = floor(ahz.o/ahz.d); HIz = LOz + ahz.n;

    /* allocate wavefield storage */
    qs = sf_complexalloc3(amx.n,amy.n,amz.n);
    qr = sf_complexalloc3(amx.n,amy.n,amz.n);
    MLOOP( qs[imz][imy][imx] = 0.0; 
	   qr[imz][imy][imx] = 0.0; );

    /* allocate image storage */
    qi = sf_floatalloc3(amx.n,amy.n,amz.n);
    MLOOP( qi[imz][imy][imx] = 0.0; );

    /* allocate cigs storage */
    qx = sf_floatalloc6(acx.n,acy.n,acz.n,ahx.n,ahy.n,ahz.n);
    HLOOP(
	CLOOP( 
	    qx[ihz-LOz][ihy-LOy][ihx-LOx][icz][icy][icx] = 0.0; 
	    );
	);
}

void img2t_init(axa amx_,
		axa amy_,
		axa amz_,
		int jcx_,
		int jcy_,
		int jcz_,
		axa aht_,
		axa aw_,
		fslice imag
    )
/*< initialize t-offset imaging condition >*/
{
    int imx,imy,imz;
    int icx,icy,icz;
    int  iht,iw;
    float ht, w;

    amx = amx_;
    amy = amy_;
    amz = amz_;

    jcx = jcx_; acx.n = amx.n / jcx;
    jcy = jcy_; acy.n = amy.n / jcy;
    jcz = jcz_; acz.n = amz.n / jcz;
    
    aht = aht_;
    aw  = aw_;

    /* allocate wavefield storage */
    qs = sf_complexalloc3(amx.n,amy.n,amz.n);
    qr = sf_complexalloc3(amx.n,amy.n,amz.n);
    MLOOP( qs[imz][imy][imx] = 0.0; 
	   qr[imz][imy][imx] = 0.0; );

    /* allocate image storage */
    qi = sf_floatalloc3(amx.n,amy.n,amz.n);
    MLOOP( qi[imz][imy][imx] = 0.0; );
    
    /* allocate cigs storage */
    qt = sf_floatalloc4(acx.n,acy.n,acz.n,aht.n);
    for( iht=0; iht<aht.n; iht++) {
	CLOOP( qt[iht][icz][icy][icx] = 0.0; );
    }

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    /* precompute phase shift */
    tt = sf_complexalloc2(aht.n,aw.n);
    for (iw=0; iw<aw.n; iw++) {
	w = aw.o+iw*aw.d;
	for (iht=0; iht<aht.n; iht++) {
	    ht = aht.o+iht*aht.d;
	    tt[iw][iht] = cexpf(-2*I*w*ht);
	}
    }
}

void img2h_init(axa amx_,
		axa amy_,
		axa amz_,
		int jcx_,
		int jcy_,
		int jcz_,
		axa ahh_,
		axa aha_,
		axa ahb_,
		axa aw_,
		fslice imag
    )
/*< initialize h-offset imaging condition >*/
{
    int imx,imy,imz;
    int icx,icy,icz;
    int  ihh;

    amx = amx_;
    amy = amy_;
    amz = amz_;

    jcx = jcx_; acx.n = amx.n / jcx;
    jcy = jcy_; acy.n = amy.n / jcy;
    jcz = jcz_; acz.n = amz.n / jcz;
    
    ahh = ahh_;
    aha = aha_;
    ahb = ahb_;
    aw  = aw_;

    /* allocate wavefield storage */
    qs = sf_complexalloc3(amx.n,amy.n,amz.n);
    qr = sf_complexalloc3(amx.n,amy.n,amz.n);
    MLOOP( qs[imz][imy][imx] = 0.0; 
	   qr[imz][imy][imx] = 0.0; );

    /* allocate image storage */
    qi = sf_floatalloc3(amx.n,amy.n,amz.n);
    MLOOP( qi[imz][imy][imx] = 0.0; );
    
    /* allocate cigs storage */
    qh = sf_floatalloc4(acx.n,acy.n,acz.n,ahh.n);
    for( ihh=0; ihh<ahh.n; ihh++) {
	CLOOP( qh[ihh][icz][icy][icx] = 0.0; );
    }
}

/*------------------------------------------------------------*/

void img2o( fslice imag,
	    fslice cigs,
	    int      iw)
/*< Apply o-offset imaging condition >*/
{
    int imx,imy,imz;
    int icx,icy,icz;

    /* imag */
    MLOOP(
	;             qi[imz][imy][imx] +=
	crealf( conjf(qs[imz][imy][imx]) 
		*     qr[imz][imy][imx]);
	);

    /* cigs */
    CLOOP(
	;             qo[icz    ][icy    ][icx    ] +=
	crealf( conjf(qs[icz*jcz][icy*jcy][icx*jcx]) 
		*     qr[icz*jcz][icy*jcy][icx*jcx] );
	);
}

void img2x( fslice imag,
	    fslice cigs,
	    int      iw)
/*< Apply x-offset imaging condition >*/
{
    int imx, imy, imz;
    int icx, icy, icz;
    int ihx, ihy, ihz;
    int imys,imyr,imzs;
    int imxs,imxr,imzr;

    /* imag */
    MLOOP(
	;             qi[imz][imy][imx] +=
	crealf( conjf(qs[imz][imy][imx]) 
		*     qr[imz][imy][imx]);
	);

    /* cigs */
/*    HLOOP(*/
/*	   XLOOP(*/
/*	       ;             qo[imz ][imy ][imx ] =*/
/*	       crealf( conjf(qs[imzs][imys][imxs]) */
/*		       *     qr[imzr][imyr][imxr]);*/
/*	       );*/
/*	   CLOOP(*/
/*	       qx[ihz-LOz][ihy-LOy][ihx-LOx][icz    ][icy    ][icx    ] +=*/
/*	       qo                           [icz*jcz][icy*jcy][icx*jcx];*/
/*	       );*/
/*	);*/

    
    for( ihz=LOz; ihz<HIz; ihz++){ 
	for( ihy=LOy; ihy<HIy; ihy++){
	    for( ihx=LOx; ihx<HIx; ihx++){

		for( icz=0; icz< acz.n; icz++){ 
		    imzs = icz*jcz - ihz,amz;
		    imzr = icz*jcz + ihz,amz;
		    if(imzs>=0 && imzs<amz.n && imzr>=0 && imzr<amz.n) {

			for( icy=0; icy< acy.n; icy++){
			    imys = icy*jcy - ihy,amy;
			    imyr = icy*jcy + ihy,amy;
			    if(imys>=0 && imys<amy.n && imyr>=0 && imyr<amy.n) {

				for( icx=0; icx< acx.n; icx++){
				    imxs = icx*jcx - ihx,amx;
				    imxr = icx*jcx + ihx,amx;
				    if(imxs>=0 && imxs<amx.n && imxr>=0 && imxr<amx.n) {

					qx[ihz-LOz][ihy-LOy][ihx-LOx][icz][icy][icx] +=
					    crealf( conjf(qs[imzs][imys][imxs]) 
						    *     qr[imzr][imyr][imxr]);

				    }
				}
			    }
			}
		    }
		}		
	    }
	}
    }
}

void img2t( fslice imag,
	    fslice cigs,
	    int      iw)
/*< Apply t-offset imaging condition >*/
{
    int imx,imy,imz,iht;
    int icx,icy,icz;
    float complex wt;

    /* imag */
    MLOOP(
	;             qi[imz][imy][imx] +=
	crealf( conjf(qs[imz][imy][imx]) 
		*     qr[imz][imy][imx] );
	);

    /* cigs */
    for(iht=0; iht<aht.n; iht++) {
	wt = tt[iw][iht];
	
	CLOOP(;             qt[iht][icz    ][icy    ][icx    ] += 
	      crealf( conjf(qs     [icz*jcz][icy*jcy][icx*jcx]) 
		      *     qr     [icz*jcz][icy*jcy][icx*jcx] * wt ); 
	    );
    }
}

void img2h( fslice imag,
	    fslice cigs,
	    int      iw)
/*< Apply h-offset imaging condition >*/
{
    int imx,imy,imz;
    int ihh,iha,ihb;
    int icx,icy,icz;

    int idx,idy,idz;
    int isx,isy,isz;
    int irx,iry,irz;

    float hh,aa,bb; /* aa,bb in radians */
    float hx,hy,hz;
    complex float cs,cr;
    float hscale;

    /* imag */
    MLOOP(
	;             qi[imz][imy][imx] +=
	crealf( conjf(qs[imz][imy][imx]) 
		*     qr[imz][imy][imx] );
	);

    /* cigs */
    for(ihh=0; ihh<ahh.n; ihh++) {                /* absolute offset */
	hh = ahh.o + ihh * ahh.d;

	hscale=1.;
	if(ihh>0) {
	    hscale=(hh+ahh.d)*aha.d;
	}

	for(ihb=0; ihb<ahb.n; ihb++) {        /* latitude  */
	    bb = ahb.o + ihb*ahb.d;
	    for(iha=0; iha<aha.n; iha++) {    /* longitude */
		aa = aha.o + iha*aha.d;
		
		hz = hh * sin(aa);
		hx = hh * cos(aa) * cos(bb);
		hy = hh * sin(aa) * sin(bb);

		/* nearest neighbour */
		idx = (int)(hx / amx.d);
		idy = (int)(hy / amy.d);
		idz = (int)(hz / amz.d);

		CLOOP(
		    isx=MM(icx*jcx-idx,amx); 
		    isy=MM(icy*jcy-idy,amy);
		    isz=MM(icz*jcz-idz,amz);
		    cs = qs[isz][isy][isx];

		    irx=MM(icx*jcx+idx,amx);
		    iry=MM(icy*jcy+idy,amy);
		    irz=MM(icz*jcz+idz,amz);
		    cr = qr[irz][iry][irx];

		    /* 
		       cs = qs @ x-hx,y-hy,z-hz
		       cr = qr @ x+hx,y+hy,z+hz
		    */

/*		    qh[ihh][icz][icy][icx] += (sqrtf(hh)/ahh.d) * crealf( conjf(cs) * cr);    */

		    qh[ihh][icz][icy][icx] += hscale * crealf( conjf(cs) * cr);    
		    ); /* cx,cy,cz */

	    } /* aa */
	}     /* bb */	
    }         /* hh */
}

/*------------------------------------------------------------*/

void img2o_close(fslice imag,
		 fslice cigs)
/*< deallocate >*/
{
    fslice_put(imag,0,qi[0][0]);
    fslice_put(cigs,0,qo[0][0]);

    img2_close();

    free(**qo); 
    free( *qo); 
    free(  qo);
}

void img2t_close(fslice imag,
		 fslice cigs)
/*< deallocate >*/
{
    fslice_put(imag,0,qi[0][0]);
    fslice_put(cigs,0,qt[0][0][0]);

    img2_close();

    free(***qt); 
    free( **qt); 
    free(  *qt); 
    free(   qt);

    free(*tt); 
    free( tt);
}

void img2h_close(fslice imag,
		 fslice cigs)
/*< deallocate >*/
{
    fslice_put(imag,0,qi[0][0]);
    fslice_put(cigs,0,qh[0][0][0]);

    img2_close();

    free(***qh); 
    free( **qh); 
    free(  *qh); 
    free(   qh);
}

void img2x_close(fslice imag,
		 fslice cigs)
/*< deallocate >*/
{
    fslice_put(imag,0,qi[0][0]);
    fslice_put(cigs,0,qx[0][0][0][0][0]);

    img2_close();

    free(*****qx); 
    free( ****qx); 
    free(  ***qx); 
    free(   **qx); 
    free(    *qx); 
    free(     qx);
}

void img2_close()
/*< deallocate >*/
{
    free(**qs); free(*qs); free(qs);
    free(**qr); free(*qr); free(qr);
    free(**qi); free(*qi); free(qi);
}
