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

#define IND(ihx,ihy) (ihy-LOy)* ahx.n + (ihx-LOx);

#define CLOOP(a) for(imy  = abs(ihy); imy<amy.n-abs(ihy) ; imy++){ \
                     imys = imy + ihy; \
                     imyr = imy - ihy; \
                 for(imx  = abs(ihx); imx<amx.n-abs(ihx) ; imx++){ \
                     imxs = imx + ihx; \
                     imxr = imx - ihx; \
                    {a} \
                 }}

static axa amx,amy,amz;
static axa ahx,ahy;
static axa aht;
static axa aw;

static float complex  **tt; /* phase shift for time offset */

static float        ****qx; /* x-offset imaging condition */
static float         ***qt; /* t-offset imaging condition */
static float          **qo; /* o-offset imaging condition */

static int LOx,HIx;
static int LOy,HIy;

/*------------------------------------------------------------*/
void imgo_init(axa amz_,
	       axa amx_,
	       axa amy_,
	       fslice imag
    )
/*< initialize o-offset imaging condition >*/
{
    int imx,imy,imz;

    amz = amz_;
    amx = amx_;
    amy = amy_;

    /* allocate image storage */
    qo = sf_floatalloc2(amx.n,amy.n);

    LOOP( qo[imy][imx] = 0.0; );
    for (imz=0; imz<amz.n; imz++) {
	fslice_put(imag,imz,qo[0]);
    }
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

    amz = amz_;
    amx = amx_;
    amy = amy_;
    aht = aht_;
    aw  = aw_;

    /* allocate image storage */
    qt = sf_floatalloc3(amx.n,amy.n,aht.n);

    for (iht=0; iht<aht.n; iht++) {
	LOOP( qt[iht][imy][imx]= 0.0; );
    }
    for (imz=0; imz<amz.n; imz++) {
	fslice_put(imag,imz,qt[0][0]);
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

void imgx_init(axa amz_,
	       axa amx_,
	       axa amy_,
	       axa ahx_,
	       axa ahy_,
	       axa aw_,
	       fslice imag
    )
/*< initialize x-offset imaging condition >*/
{
    int imx,imy,imz;
    int ihx,ihy;

    amz = amz_;
    amx = amx_;
    amy = amy_;
    ahx = ahx_;
    ahy = ahy_;
    aw  = aw_;

    /* allocate image storage */
    qx = sf_floatalloc4(amx.n,amy.n,ahx.n,ahy.n);

    for (ihy=0; ihy<ahy.n; ihy++) {
	for (ihx=0; ihx<ahx.n; ihx++) {
	    LOOP( qx[ihy][ihx][imy][imx]= 0.0; );
	}
    }
    for (imz=0; imz<amz.n; imz++) {
	fslice_put(imag,imz,qx[0][0][0]);
    }

    LOx = floor(ahx.o/ahx.d); HIx = LOx + ahx.n;
    LOy = floor(ahy.o/ahy.d); HIy = LOy + ahy.n;

}

/*------------------------------------------------------------*/

void imgo(int            imz,
	  int            iw,
	  fslice         imag  /* image slice [ny][nx] */,
	  float complex **ww_s /* source   wavefield */,
	  float complex **ww_r /* receiver wavefield */
    )
/*< Apply o-offset imaging condition >*/
{
    int imx,imy;

    fslice_get(imag,imz,qo[0]);
    LOOP(;              qo  [imy][imx] += 
	  crealf( conjf(ww_s[imy][imx]) 
		  *     ww_r[imy][imx] ); );
    fslice_put(imag,imz,qo[0]);
}

void imgt(int            imz,
	  int            iw,
	  fslice         imag  /* image slice [ny][nx] */,
	  float complex **ww_s /* source   wavefield */,
	  float complex **ww_r /* receiver wavefield */
    )
/*< Apply t-offset imaging condition >*/
{
    int imx,imy,iht;

    fslice_get(imag,imz,qt[0][0]);
    for(iht=0; iht<aht.n; iht++) {
	LOOP(;          qt[iht][imy][imx] += 
	     crealf( conjf(ww_s[imy][imx]) 
		     *     ww_r[imy][imx] * tt[iw][iht] ); );
    }
    fslice_put(imag,imz,qt[0][0]);
}

void imgx(int            imz,
	  int            iw,
	  fslice         imag  /* image slice [ny][nx] */,
	  float complex **ww_s /* source   wavefield */,
	  float complex **ww_r /* receiver wavefield */
    )
/*< Apply x-offset imaging condition >*/
{
    int imx,imy,ihx,ihy;
    int imys,imyr;
    int imxs,imxr;

    fslice_get(imag,imz,qx[0][0][0]);

    for(ihy=LOy; ihy<HIy; ihy++){
	for(ihx=LOx; ihx<HIx; ihx++){
	    CLOOP(
		qx[ihy-LOy][ihx-LOx][imy][imx] += crealf( 
		    conjf(ww_s[imys][imxs]) 
		    *     ww_r[imyr][imxr]          ); 
		);
	}
    }
    fslice_put(imag,imz,qx[0][0][0]);
}

/*------------------------------------------------------------*/

void imgo_close()
/*< deallocate >*/
{
    free( *qo); free( qo);
}

void imgt_close()
/*< deallocate >*/
{
    free( *tt); free( tt);
    free(**qt); free( *qt); free( qt);
}

void imgx_close()
/*< deallocate >*/
{
    free(***qx); free(**qx); free( *qx); free( qx);
}
