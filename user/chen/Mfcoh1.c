/* Fast C1 Coherence */

/*
  Copyright (C) 2013 Zhonghuan Chen, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranyw of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>
#include "comp.h"
#include "polyfit.h"


static void *h1, *h2;

int main(int argc, char* argv[])
{
	sf_file in, out, idip, xdip;
	int n1, n2, n3, ntw, m1;
	int i1, i2, i3, j1, k1, k2, ix, iy, it;
	float ***u1, ***u2, **u3,  *c1, *c2;
	float max, ***v1, ***v2, d1, d2, c[3];
	bool twod, verb;
	int min1, min2, max1, max2;

	sf_init(argc, argv);

	in  = sf_input("in");	/* 3D data set */
	out = sf_output("out");	/* 3D coherence cube */

	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histint(in,"n3",&n3)) n3=1; 

	if (!sf_getint("ntw", &ntw) ) ntw=5;
	/* half window size for coherence */
	if (!sf_getbool("twod",&twod) ) twod=true;
	/* y: only twod coherence */
	if (!sf_getbool("verb",&verb) ) verb=true;
	/* verbosity */
	if(sf_getstring("idip")!=NULL) idip  = sf_output("idip");
	/* inline dip */
	else idip = NULL;
	if(sf_getstring("xdip")!=NULL) xdip  = sf_output("xdip");
	/* crossline dip */
	else xdip = NULL;
	if (!sf_getint("min1", &min1) ) min1=-2;
	if (!sf_getint("max1", &max1) ) max1=2;
	/* inline slope */
	if (!sf_getint("min2", &min2) ) min2=-2;
	if (!sf_getint("max2", &max2) ) max2=2;
	/* xline slope */


	h1 = polyfit_init(max1-min1+1, 3, min1, 0);

	if(n3==1)  twod=true;

	u1 = sf_floatalloc3(n1, n2, n3);
	u2 = sf_floatalloc3(n1, n2, n3);
	if(idip)v1 = sf_floatalloc3(n1, n2, n3);

	m1 = 2*ntw+1;
	u3 = sf_floatalloc2(m1, 2);
	c1 = sf_floatalloc(max1-min1+1);

	if(twod) {min2=0; max=2; h2=NULL;}
	else {
		h2 = polyfit_init(max2-min2+1, 3, min2, 0);
		if(xdip)v2 = sf_floatalloc3(n1, n2, n3);
		c2 = sf_floatalloc(max2-min2+1);
	}

	sf_floatread(u1[0][0], n1*n2*n3, in);
	for(i3=0; i3<n3; i3++)
	{
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++)
		{
			max = 0.0;
			for(j1=-ntw; j1<=ntw; j1++)
			{
				it =i1 + j1;
				if(it<0 || it>n1-1) u3[0][j1+ntw] = 0.0;
				else u3[0][j1+ntw] = u1[i3][i2][it];
			}

			for(k1=min1; k1<=max1; k1++)
			{
				ix = i2 + 1;
				if(ix>= n2 || ix <0) c1[k1-min1] = 0.0;
				else {
					for(j1=-ntw; j1<=ntw; j1++)
					{
						it = i1 + k1 + j1;
						if(it <0 || it>=n1)	u3[1][j1+ntw] = 0.0;
						else u3[1][j1+ntw] = u1[i3][ix][it];
					}
					c1[k1-min1] = comp_ncc(u3[0], u3[1], m1);
				}
			}
			polyfit_coef(h1, c1, c);
			d1 = -0.5*c[1]/c[2];
			max = c[0] + 0.5 * d1*c[1];
			u2[i3][i2][i1] = max;
			if(idip) v1[i3][i2][i1] = d1;
			
			if(!twod){
				for(k2=min2; k2<=max2; k2++)
				{
					iy = i3 + 1;
					if(iy>= n3 || iy <0) c2[k2-min2] = 0.0;
					else {
						for(j1=-ntw; j1<=ntw; j1++)
						{
							it = i1 + k2 + j1;
							if(it <0 || it>=n1)	u3[1][j1+ntw] = 0.0;
							else u3[1][j1+ntw] = u1[iy][i2][it];
						}	
						c2[k2-min2] = comp_ncc(u3[0], u3[1], m1);
					}
				}
				polyfit_coef(h2, c2, c);
				d2 = -0.5*c[1]/c[2];
				max = c[0] + 0.5 * d2*c[1];
				u2[i3][i2][i1] *= max;
				if(xdip) v2[i3][i2][i1] = d2;
			}
		}
		if(verb) sf_warning("%d of %d;", i3, n3);
	}
	sf_floatwrite(u2[0][0], n1*n2*n3, out);
	if(idip)sf_floatwrite(v1[0][0], n1*n2*n3, idip);
	if(!twod)
	{
		free(c2);
		polyfit_close(h2);
		if(xdip)
		{
			sf_floatwrite(v2[0][0], n1*n2*n3, xdip);
			free(**v2);
			free(*v2);
			free(v2);
		}
	}
	free(c1);
	polyfit_close(h1);
	free(**u1);
	free(*u1);
	free(u1);
	free(**u2);
	free(*u2);
	free(u2);
	free(*u3);
	free(u3);
	if(idip)
	{
		free(**v1);
		free(*v1);
		free(v1);
	}
	return 0;
}






