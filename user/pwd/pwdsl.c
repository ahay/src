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
#include "pwsmooth.h"
#include "repeat2.h"

static int n;
static float *tmp;

void pwdsl_init(int n1, int n2       /* data dimensions */,
		int order            /* accuracy order */,
		int rect1, int rect2 /* triangle radius */,
		float eps            /* regularization parameter */)                
/*< initialize >*/
{
    pwsmooth_init (rect2,n1,n2,order,eps);
    sf_triangle1_init(rect1,n1);

    repeat2_init(n1,n2,sf_triangle1_lop);

    n = n1*n2;
    tmp = sf_floatalloc(n);
}

void pwdsl_set(float **dip /* local slope */)
/*< set local slope >*/
{
    pwsmooth_set(dip);
}

void pwdsl_close(void)
/*< free allocated storage >*/
{
    sf_triangle1_close();
    pwsmooth_close();
    
    free(tmp);
}

void pwdsl_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    if (n != nx || n != ny) sf_error("__FILE__: wrong size");
    
    sf_chain(pwsmooth_lop,repeat2_lop,adj,add,n,n,n,x,y,tmp);    
}
