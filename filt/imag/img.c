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

static axa amx,amy,amz;

/*------------------------------------------------------------*/

void img_init(axa amz_        /* depth */,
	      axa amx_        /* i-line (data) */,
	      axa amy_        /* x-line (data) */
    )
/*< initialize SR imaging condition >*/
{
    amz = amz_;
    amx = amx_;
    amy = amy_;
}

/*------------------------------------------------------------*/

void img_zero(fslice imag,
	      float         **qq
    )
/*< Zero output image file >*/
{
    int imx,imy,imz;

    LOOP( qq[imy][imx] = 0.0; );
    for (imz=0; imz<amz.n; imz++) {
	fslice_put(imag,imz,qq[0]);
    }
}

/*------------------------------------------------------------*/

void img_xcor(int            iz,
	      fslice         imag  /* image slice [nz][ny][nx] */,
	      float         **qq   /* image array */,
	      float complex **ww_s /* source   wavefield */,
	      float complex **ww_r /* receiver wavefield */
    )
/*< Apply x-correlation imaging condition >*/
{
    int imx,imy;

    fslice_get(imag,iz,qq[0]); /* imaging */
    LOOP(;             qq  [imy][imx] += 
	 crealf( conjf(ww_s[imy][imx]) * ww_r[imy][imx] ); );
    fslice_put(imag,iz,qq[0]);
}

/*------------------------------------------------------------*/

