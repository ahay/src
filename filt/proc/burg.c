/* PEF estimation by Burg's algorithm. */
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

#include "burg.h"

float pef_burg(int n, float* trace /* data [n] */)
/*< return a for 3-point pef (1,a,1) >*/
{
    int i;
    float a, avto, cros;
    
    avto = trace[n-1]*trace[n-1];
    a = avto + trace[0]*trace[0];
    cros = 0.;
    for (i=0; i < n-1; i++) {
	avto += trace[i]*trace[i];
	cros += trace[i]*trace[i+1];
    }
    if (avto == 0.) { 
	a=-2.;
    } else {
	a = - (a + 2.*cros)/avto;
    }

    return a;
}

float pef_burg2(int n1, int n2, float** trace /* data[n2][n1] */)
/*< return a for 3-point pef (1,a,1) >*/
{
    int i1, i2;
    float a;
    double avto, cros;

    a = avto = cros = 0.;
    
    for (i2=0; i2 < n2; i2++) {
	avto += trace[i2][n1-1]*trace[i2][n1-1];
	a += (trace[i2][n1-1]*trace[i2][n1-1] + trace[i2][0]*trace[i2][0]);
	for (i1=0; i1 < n1-1; i1++) {
	    avto += trace[i2][i1]*trace[i2][i1];
	    cros += trace[i2][i1]*trace[i2][i1+1];
	}
    }

    if (avto == 0.) { 
	a=-2.;
    } else {
	a = - (a + 2.*cros)/avto;
    }
    
    return a;
}

void pef_define (int n, float a, float eps, float* diag, float** offd)
/*< insert PEF into a pentadiagonal matrix for regularization >*/
{
    int i;

    for (i=0; i < n-2; i++) {
	diag[i+1] = (2. + a*a)*eps; 
	offd[0][i] = (2.*a)*eps;
	offd[1][i] = eps;
    }
    offd[0][0]   = (1.+2.*a)*eps;
    offd[0][n-2] = (1.+2.*a)*eps;
    diag[0]    = (1.+(1.+a)*(1.+a))*eps;
    diag[n-1]  = (1.+(1.+a)*(1.+a))*eps;
}
