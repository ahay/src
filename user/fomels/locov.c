/* Local covariance filter */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "locov.h"

static float f[3];

void locov_init(float a, float b)
/*< initialize >*/
{
    float ab, ab1;

    ab = a*b;
    ab1 = (a-1.)*b;

    f[0] = 1+ab*(2.5+3.*ab1);
    f[1] = -ab*(4./3.+2.*ab1);
    f[2] = ab*(1./12+0.5*ab1);
}

void locov_filt(int n, const float *inp, float *out)
/*< filter >*/
{
    int i;

    out[0] = f[0]*inp[0] + f[1]*(inp[0]+inp[1]) + f[2]*(inp[0]+inp[2]);
    out[1] = f[0]*inp[1] + f[1]*(inp[0]+inp[2]) + f[2]*(inp[0]+inp[3]);
    for (i=2; i < n-2; i++) {
	out[i] = f[0]*inp[i] + f[1]*(inp[i-1]+inp[i+1]) + f[2]*(inp[i-2]+inp[i+2]);
    }
    out[n-2] = f[0]*inp[n-2] + f[1]*(inp[n-3]+inp[n-1]) + f[2]*(inp[n-4]+inp[n-1]);
    out[n-1] = f[0]*inp[n-1] + f[1]*(inp[n-2]+inp[n-1]) + f[2]*(inp[n-3]+inp[n-1]);
}

