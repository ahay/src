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

#include "img.h"

#include "slice.h"
/*^*/

#define LOOP(a) for( imy=0; imy< amy.n; imy++){ \
                for( imx=0; imx< amx.n; imx++){ {a} }}

#define MLOOP(a) for( imz=0; imz< amz.n; imz++){ \
                 for( imy=0; imy< amy.n; imy++){ \
                 for( imx=0; imx< amx.n; imx++){ {a} }}}

#define HLOOP(a) for(ihz=LOz; ihz<HIz; ihz++){ \
                 for(ihy=LOy; ihy<HIy; ihy++){ \
                 for(ihx=LOx; ihx<HIx; ihx++){ \
                    {a} \
                 }}}

#define CLOOP(a) for(imz  = abs(ihz); imz<amz.n-abs(ihz) ; imz++){ \
                     imzs = imz + ihz; \
                     imzr = imz - ihz; \
                 for(imy  = abs(ihy); imy<amy.n-abs(ihy) ; imy++){ \
                     imys = imy + ihy; \
                     imyr = imy - ihy; \
                 for(imx  = abs(ihx); imx<amx.n-abs(ihx) ; imx++){ \
                     imxs = imx + ihx; \
                     imxr = imx - ihx; \
                    {a} \
                 }}}

#define IND(ihx,ihy,ihz) (ihz-LOz)*(ahx.n*ahy.n) + \
                         (ihy-LOy)* ahx.n        + \
                         (ihx-LOx);

static axa amx,amy,amz;
static axa ahx,ahy,ahz;
static axa aht;
static axa aw;

static float complex  **tt; /* phase shift for time offset */
static float complex ***qs,***qr;
static float         ***qi;

static int LOx,HIx;
static int LOy,HIy;
static int LOz,HIz;

/*------------------------------------------------------------*/

void imgstore( int imz,
	       float complex **ww_s,
	       float complex **ww_r
    )
/*< store wavefield >*/
{
    int imx,imy;

    LOOP( 
	qs[imz][imy][imx] = ww_s[imy][imx];
	qr[imz][imy][imx] = ww_r[imy][imx];
	);
}

/*------------------------------------------------------------*/
void imgo_init(axa amz_,
	       axa amx_,
	       axa amy_,
	       fslice imag
    )
/*< initialize o-offset imaging condition >*/
{
    int imx,imy,imz;

    amx = amx_;
    amy = amy_;
    amz = amz_;

    /* allocate image storage */
    qs = sf_complexalloc3(amx.n,amy.n,amz.n);
    qr = sf_complexalloc3(amx.n,amy.n,amz.n);
    qi = sf_floatalloc3  (amx.n,amy.n,amz.n);
    MLOOP( 
	qs[imz][imy][imx] = 0.0; 
	qr[imz][imy][imx] = 0.0; 
	qi[imz][imy][imx] = 0.0; 
	);
    fslice_put(imag,1,qi[0][0]);
}

void imgx_init(axa amz_,
	       axa amx_,
	       axa amy_,
	       axa ahx_,
	       axa ahy_,
	       axa ahz_,
	       fslice imag
    )
/*< initialize x-offset imaging condition >*/
{
    int imx,imy,imz;
    int ihx,ihy,ihz,ih;

    amx = amx_;
    amy = amy_;
    amz = amz_;

    ahx = ahx_;
    ahy = ahy_;
    ahz = ahz_;

    LOx = floor(ahx.o/ahx.d); HIx = LOx + ahx.n;
    LOy = floor(ahy.o/ahy.d); HIy = LOy + ahy.n;
    LOz = floor(ahz.o/ahz.d); HIz = LOz + ahz.n;

    /* allocate image storage */
    qs = sf_complexalloc3(amx.n,amy.n,amz.n);
    qr = sf_complexalloc3(amx.n,amy.n,amz.n);
    qi = sf_floatalloc3  (amx.n,amy.n,amz.n);
    MLOOP( 
	qs[imz][imy][imx] = 0.0; 
	qr[imz][imy][imx] = 0.0; 
	qi[imz][imy][imx] = 0.0; 
	);

    HLOOP( ih = IND(ihx,ihy,ihz);
	   fslice_put(imag,ih,qi[0][0]);
	);
}

void imgt_init(axa amz_,
	       axa amx_,
	       axa amy_,
	       axa aht_,
	       axa aw_,
	       fslice imag
    )
/*< initialize t-offset imaging condition >*/
{
    int imx,imy,imz;
    int  iht,iw;
    float ht, w;

    amx = amx_;
    amy = amy_;
    amz = amz_;

    aht = aht_;
    aw  = aw_;

    /* allocate image storage */
    qs = sf_complexalloc3(amx.n,amy.n,amz.n);
    qr = sf_complexalloc3(amx.n,amy.n,amz.n);
    qi = sf_floatalloc3  (amx.n,amy.n,amz.n);
    MLOOP( 
	qs[imz][imy][imx] = 0.0; 
	qr[imz][imy][imx] = 0.0; 
	qi[imz][imy][imx] = 0.0; 
	);

    for (iht=0; iht<aht.n; iht++) {
	fslice_put(imag,iht,qi[0][0]);
    }
    
    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    /* phase shift */
    tt = sf_complexalloc2(aht.n,aw.n);
    for (iw=0; iw<aw.n; iw++) {
	w = aw.o+iw*aw.d;
	for (iht=0; iht<aht.n; iht++) {
	    ht = aht.o+iht*aht.d;
	    tt[iw][iht] = cexpf(-2*I*w*ht);
	}
    }
}
/*------------------------------------------------------------*/

void imgo( fslice imag,
	      int   iw)
/*< Apply o-offset imaging condition >*/
{
    int imx,imy,imz;

    fslice_get(imag,0,qi[0][0]);
    MLOOP(
	;             qi[imz][imy][imx] +=
	crealf( conjf(qs[imz][imy][imx]) 
		*     qr[imz][imy][imx] );
	);
    fslice_put(imag,0,qi[0][0]);
}

void imgx( fslice imag,
	      int   iw)
/*< Apply x-offset imaging condition >*/
{
    int imx ,imy ,imz;
    int ihx ,ihy ,ihz ,ih;
    int imys,imyr,imzs;
    int imxs,imxr,imzr;

    HLOOP( ih = IND(ihx,ihy,ihz);
	   fslice_get(imag,ih,qi[0][0]);
	   CLOOP(
	       ;             qi[imz ][imy ][imx ] +=
	       crealf( conjf(qs[imzs][imys][imxs]) 
		       *     qr[imzr][imyr][imxr] );
	       );
	   fslice_put(imag,ih,qi[0][0]);
	);
}

void imgt( fslice imag,
	   int      iw)
/*< Apply t-offset imaging condition >*/
{
    int imx,imy,imz,iht;
    float complex wt;

    for(iht=0; iht<aht.n; iht++) {
	wt = tt[iw][iht];

	fslice_get(imag,iht,qi[0][0]);
	MLOOP(;             qi[imz][imy][imx] += 
	      crealf( conjf(qs[imz][imy][imx]) 
		      *     qr[imz][imy][imx] * wt ); );
	fslice_put(imag,iht,qi[0][0]);
    }
}

/*------------------------------------------------------------*/

void imgo_close()
/*< deallocate >*/
{
}

void imgt_close()
/*< deallocate >*/
{
    free( *tt); free( tt);
}

void imgx_close()
/*< deallocate >*/
{
}
