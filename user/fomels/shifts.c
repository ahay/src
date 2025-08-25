/*
  Copyright (C) 2025 University of Texas at Austin
  
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

void shifts2(int s1, int s2, int n1, int n2, float **inp, float ***sft)
/*< generate shifts in 2D >*/
{
    int i1, i2, j1, j2, i, k1, k2, l1, l2;
    float a;
    
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    i=0;
	    for (j2=0; j2 <= s2; j2++) {
		for (j1=-s1; j1 <= s1; j1++) {
		    a = 0.0f;
		    k2 = i2+j2;
		    k1 = i1+j1;
		    if (k1 >=0 && k1 < n1 && k2 >=0 && k2 < n2) {
			a += inp[k2][k1];
		    }
		    l2 = i2-j2;
		    l1 = i1-j1;
		    if (l1 >=0 && l1 < n1 && l2 >=0 && l2 < n2) {
			a += inp[l2][l1];
		    }
		    sft[i][i2][i1] = a;
		    i++;
		}
	    }
	}
    }
}

void shifts1(int s1, int n1, float *inp, float **sft)
/*< generate shifts in 1D >*/
{
    int i1, j1, k1;
    
    for (i1=0; i1 < n1; i1++) {
	for (j1=-s1; j1 <= s1; j1++) {
	    k1 = i1+j1;
	    if (k1 >=0 && k1 < n1) {
		sft[j1+s1][i1] = inp[k1];
	    }
	}
    }
}
