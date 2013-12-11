/* 1-D missing data interpolation with known filter */
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

#include <rsf.h>
/*^*/

#include "nfmis.h"

#include "nfconv.h"

//static int nx, ny;
//static float *zero;



//////////////////////////////////////////////////////////////////

void nfmis(int niter         /* number of iterations */,
      int n1, int na2, 
      sf_complex *filt,
	  sf_complex *xx         /* data/model */, 
	  const bool *known /* mask for known data */,
	  bool verb         /* verbosity flag */) 
/*< interpolate >*/
{    
    int ix;
    sf_complex *zero;
    zero = sf_complexalloc(n1);
    for (ix=0; ix < n1; ix++) {
	zero[ix]=sf_cmplx(0,0);
    }

    nfconv_init(filt, n1, na2);

    //sf_warning("xxxx[0]=%f",crealf(xx[0]));
    //sf_ccdstep_init();
	sf_csolver (nfconv_lop, sf_ccgstep, n1, n1, xx, zero, 
		       niter, "x0", xx, "known", known, "verb",verb,"end");


    free(zero);
	sf_ccgstep_close();




}

////////////////////////////////////////////////////////////////////////





/* 	$Id: mis1.c 7107 2011-04-10 02:04:14Z ivlad $	 */
