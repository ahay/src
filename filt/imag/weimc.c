/* 3-D imaging conditions for shot-profile WE migration */
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

#include "weimc.h"

#include "slice.h"
/*^*/


#define MLOOP(a) for(imz=0; imz<amz.n; imz++){ \
                 for(imy=0; imy<amy.n; imy++){ \
                 for(imx=0; imx<amx.n; imx++){ \
                    {a} \
                 }}}

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

static bool verb;
static axa amx,amy,amz, aw;
static axa ahx,ahy,ahz;

static float         ***qq;
static float complex ***us;
static float complex ***ur;

static int LOx,HIx;
static int LOy,HIy;
static int LOz,HIz;

/*------------------------------------------------------------*/

void weimc_init(bool verb_,
		axa amx_        /* i-line (data) */,
		axa amy_        /* x-line (data) */,
		axa amz_        /* depth */,
		axa aw_         /* frequency */
)
/*< initialize WE IC >*/
{
    verb=verb_;

    amz = amz_;
    amx = amx_;
    amy = amy_;
    aw  = aw_;

    /* allocate wavefield storage */
    us = sf_complexalloc3(amx.n,amy.n,amz.n);
    ur = sf_complexalloc3(amx.n,amy.n,amz.n);
    qq = sf_floatalloc3  (amx.n,amy.n,amz.n);
}

void hhimc_init(axa ahx_,
		axa ahy_,
		axa ahz_
)
/*< initialize prestack IC >*/
{
    ahx = ahx_;
    ahy = ahy_;
    ahz = ahz_;

    LOx = floor(ahx.o/ahx.d); HIx = LOx + ahx.n;
    LOy = floor(ahy.o/ahy.d); HIy = LOy + ahy.n;
    LOz = floor(ahz.o/ahz.d); HIz = LOz + ahz.n;
    
    sf_warning("LOx=%d HIx=%d",LOx,HIx);
    sf_warning("LOy=%d HIy=%d",LOy,HIy);
    sf_warning("LOz=%d HIz=%d",LOz,HIz);

}
/*------------------------------------------------------------*/

void weimc_close(void)
/*< free allocated storage >*/
{
    free(**us); free( *us); free( us);
    free(**ur); free( *ur); free( ur);
    free(**qq); free( *qq); free( qq);
}

/*------------------------------------------------------------*/

void zoimc(fslice sdat /* source   data [nw][nz][ny][nx] */,
	   fslice rdat /* receiver data [nw][nz][ny][nx] */,
	   fslice imag /*         image     [nz][ny][nx] */
    )
/*< apply imaging condition >*/
{
    int imx,imy,imz,iw;

    MLOOP( qq[imz][imy][imx] = 0.0; );
    fslice_put(imag,0,qq[0][0]);
	
    for (iw=0; iw<aw.n; iw++) {
	if(verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);
	
	fslice_get(sdat,iw,us[0][0]);
	fslice_get(rdat,iw,ur[0][0]);
	fslice_get(imag, 0,qq[0][0]);	    
	MLOOP(;         qq[imz][imy][imx] += crealf( 
		  conjf(us[imz][imy][imx]) 
		  *     ur[imz][imy][imx]          ); 
	    );
	fslice_put(imag, 0,qq[0][0]); 
    } /* w */
}

/*------------------------------------------------------------*/

void hhimc(fslice sdat /* source   data [nw][nz][ny][nx] */,
	   fslice rdat /* receiver data [nw][nz][ny][nx] */,
	   fslice imag /*         image     [nz][ny][nx] */
    )
/*< apply imaging condition >*/
{
    int iw,ih;
    int imx, imy, imz;
    int ihx, ihy, ihz;
    int imxs,imys,imzs;
    int imxr,imyr,imzr;

    MLOOP( qq[imz][imy][imx] = 0.0; );
    HLOOP( ih = IND(ihx,ihy,ihz);
	   fslice_put(imag,ih,qq[0][0]);
	);
    
    for (iw=0; iw<aw.n; iw++) {
	if(verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);
	
	fslice_get(sdat,iw,us[0][0]);
	fslice_get(rdat,iw,ur[0][0]);

	HLOOP( ih = IND(ihx,ihy,ihz);
	       fslice_get(imag,ih,qq[0][0]);
	       
	       CLOOP(;         qq[imz ][imy ][imx ] += crealf( 
			 conjf(us[imzs][imys][imxs]) 
			 *     ur[imzr][imyr][imxr]          ); 
		   );

	       fslice_put(imag,ih,qq[0][0]);
	    );
    } /* w */
}

