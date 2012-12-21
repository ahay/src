/* Gaussian derivative models */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>
#include "hermite.h"

static int n;
static float sigma;

void dgauss_init(int m, float s)
/*< initialize >*/
{
	n = m;
	sigma = fabs(s);
	hermite_init(m);
}

void dgauss_close()
/*< release memory >*/
{
	hermite_close();
}


float dgauss_val(int k, float x)
/*< Polynomial value >*/
{
	float r1, r2;
	r2 = x/sigma;
	r1 = hermite_val(k, r2);
	r2 = r1*exp(-r2*r2/2.0)/sqrt(2.0*SF_PI)/sigma;
	if (k%2 == 1) return -r2;
	else return r2;
}



