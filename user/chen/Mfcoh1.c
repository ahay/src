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
#include "fcoh1.h"

int main(int argc, char* argv[])
{
	sf_file in, out, idip, xdip;
	int n1, n2, n3, ntw;
	int i1, i2, i3;
	float **u1, **u2, **pu, **p1, **p2;
	bool twod, verb;
	int min1, min2, max1, max2;
	void *h1, *h2;

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


	if(n3==1)  twod=true;

	u1 = sf_floatalloc2(n1, n2);
	p1 = sf_floatalloc2(n1, n2);

	h1 = fcoh1_init(n1, min1, max1, ntw);

	if(twod) for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);

		/* initial for first trace */
		fcoh1_acorr(u1[0], p1[0], n1, ntw);

		for(i2=0; i2<n2-1; i2++)
		{
			fcoh1_acorr(u1[i2+1], p1[i2+1], n1, ntw);
			fcoh1_tr(h1, u1[i2], u1[i2+1], p1[i2], p1[i2+1], NULL);
		}
		for(i1=0; i1<n1; i1++)
		{
			u1[n2-1][i1] = u1[n2-2][i1];
			if (idip) p1[n2-1][i1] = p1[n2-2][i1];
		}
		sf_floatwrite(u1[0], n1*n2, out);
		if(idip) sf_floatwrite(p1[0], n1*n2, idip);
	}else {  /* end 2D begin 3D */
		u2 = sf_floatalloc2(n1, n2);
		p2 = sf_floatalloc2(n1, n2);

		h2 = fcoh1_init(n1, min2, max2, ntw);

		sf_floatread(u2[0], n1*n2, in);
		for(i2=0; i2<n2; i2++)
			fcoh1_acorr(u2[i2], p2[i2], n1, ntw);
		for(i3=1; i3<n3; i3++)
		{
			sf_floatread(u1[0], n1*n2, in);

			fcoh1_acorr(u1[0], p1[0], n1, ntw);
			fcoh1_tr(h2, u2[0], u1[0], p2[0], p1[0], NULL);

			for(i2=1; i2<n2; i2++)
			{
				fcoh1_acorr(u1[i2], p1[i2], n1, ntw);
				fcoh1_tr(h2, u2[i2], u1[i2], p2[i2], p1[i2], NULL);
				fcoh1_tr(h1, u1[i2-1], u1[i2], p1[i2-1], p1[i2], u2[i2-1]);
			}
			for(i1=0; i1<n1; i1++)
			{
				u2[n2-1][i1] = u2[n2-1][i1];
				p1[n2-1][i1] = p1[n2-2][i1];
			}
			pu = u2; u2 = u1; u1 = pu;
			sf_floatwrite(u2[0], n1*n2, out);
			if(idip) sf_floatwrite(p1[0], n1*n2, idip);
			if(xdip) sf_floatwrite(p2[0], n1*n2, xdip);
			if(verb) sf_warning("%d of %d;", i3, n3);
		}
		sf_floatwrite(u2[0], n1*n2, out);
		free(*u2);
		free(u2);
		free(*p2);
		free(p2);
		fcoh1_close(h2);
	}
	fcoh1_close(h1);
	free(*u1);
	free(u1);
	free(*p1);
	free(p1);
	return 0;
}






