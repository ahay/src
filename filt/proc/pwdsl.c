/* Triangle smoothing along slope using plane-wave destruction */
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

#include "pwdsl.h"
#include "predict.h"
#include "repeat.h"

static float *tmp1, *tmp2;

void pwdsl_init(int m1, int m2       /* data dimensions */, 
		int rect1, int rect2 /* triangle radius */,
		float eps            /* regularization parameter */)
/*< initialize >*/
{
    predict_init (m1,m2,eps,rect2);
    sf_triangle1_init(rect1,m1);
    repeat_init(m1,m2,sf_triangle1_lop);

    tmp1 = sf_floatalloc(m1*m2);
    tmp2 = sf_floatalloc(m1*m2);
}

void pwdsl_set(float **dip /* dip field [n2][n1] */)
/*< set the local slopes for applying the linear operator >*/
{
    predict_set(dip);
}

void pwdsl_close(void)
/*< free allocated storage >*/
{
    sf_triangle1_close();
    predict_close();
    
    free(tmp1);
    free(tmp2);
}

void pwdsl_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    float *inp, *out;

    if (adj) {
	inp = y;
	out = x;
    } else {
	inp = x;
	out = y;
    }

    predict_lop  (true,  false, nx, nx, tmp1, inp);
    subtract_lop (true,  false, nx, nx, tmp2, tmp1);
    repeat_lop   (true,  false, nx, nx, tmp1, tmp2);
    repeat_lop   (false, false, nx, nx, tmp1, tmp2);
    subtract_lop (false, false, nx, nx, tmp2, tmp1);
    predict_lop  (false, add,   nx, nx, tmp1, out);
}
