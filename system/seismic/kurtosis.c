/* Compute kurtosis */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

float kurtosis(int n, float *data /* [n] */)
/*< output kurtosis >*/
{
    int i;
    float den, num, kur, s2, s4;

    num = den = 0.;
    for (i=0; i < n; i++) {
	s2 = data[i];
	s2 = s2*s2;
	s4 = s2*s2;
	num += s4;
	den += s2;
    }
    kur = n*num/(den*den);

    return kur;
}
