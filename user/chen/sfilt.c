/* simple filters  */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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

typedef float (*sfilt)(int,float*);
/* generic filter interface */
/*^*/


float sfilt_mean(int n, float *x)
/*< mean filter >*/
{
	float d=0.0;
	do{
		n--;
		d += x[n-1];
	}while(n>0);
	return (d/n);
}

float sfilt_median(int n, float *p)
/*< median filter >*/
{
	int i1, j1, chg;
	float temp;

	for(j1=n-1; j1>=n/2; j1--)
	{
		chg=0;
		for(i1=0; i1<j1; i1++)
		if(p[i1] > p[i1+1])
		{
			temp = p[i1];
			p[i1] = p[i1+1];
			p[i1+1] = temp;
			chg=1;
		}
		if(chg==0) break;
	}
	return (p[n/2]);
}



