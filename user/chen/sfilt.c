/* simple filters  */

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

float sfilt_mean(int n, float *x)
/*< mean filter >*/
{
	int i;
	double d=0.0;
	for (i=0; i<n; i++) d += x[i];
	return (d/n);
}

float sfilt_median(int n, float *p)
/*< median filter >*/
{
	int i1, j1, chg;
	float temp;

	for(j1=n-1; j1>=n/2+1; j1--)
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


typedef float (*sfilt)(int,float*);
/* generic filter interface */
/*^*/

sfilt sfilt_c2f(char *c)
/*< filter selecter >*/
{
	if(strcmp(c, "mean") == 0) return sfilt_mean;
	else if(strcmp(c, "median") == 0) return sfilt_median;
	else return NULL;
}


