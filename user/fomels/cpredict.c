/* Trace prediction with plane-wave destruction for constant dips */
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

#include "cpredict.h"
#include "cpwd.h"

static int n1, nb;
static sf_bands slv;
static float *offd, eps, *tt;
static cpwd w;

void cpredict_init (int nx  /* data size */,
		    float e /* regularization parameter */)
/*< initialize >*/
{
    const int nw=1;

    n1 = nx;
    nb = 2*nw;

    eps = e;

    slv = sf_banded_init (n1, nb);
    offd = sf_floatalloc (nb);

    w = cpwd_init (n1, nw);
    tt = sf_floatalloc (n1);
}

void cpredict_close (void)
/*< free allocated storage >*/
{
    sf_banded_close (slv);
    free (offd);
    cpwd_close (w);
    free (tt);
}

void cpredict_step(bool adj     /* adjoint flag */,
		   bool forw    /* forward or backward */, 
		   float* trace /* trace */, 
		   float pp    /* slope */)
/*< prediction step >*/
{
    float diag;
    
    offd[0] = -4.*eps;
    offd[1] = eps;
    diag = 6.*eps + cpwd_define (forw, w, pp, offd);
    sf_banded_const_define (slv, diag, offd);

    if (adj) {
	sf_banded_solve (slv, trace);
	cpwd_set (true, w, trace, trace, tt);
    } else {
	cpwd_set (false, w, trace, trace, tt);
	sf_banded_solve (slv, trace);
    }
}
