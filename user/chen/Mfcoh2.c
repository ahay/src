/* Fast C2 Coherence */

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
#include "fcoh2.h"

int main(int argc, char* argv[])
{
	sf_file in, out, idip;
	int n1, n2, n3, ntw, nxw;
	int i3;
	float **u1, **p1;
	bool twod, verb;
	int min1, min2, max1, max2;
	void *h1;

	sf_init(argc, argv);

	in  = sf_input("in");	/* 3D data set */
	out = sf_output("out");	/* 3D coherence cube */

	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histint(in,"n3",&n3)) n3=1; 

	if (!sf_getint("ntw", &ntw) ) ntw=5;
	/* half window size for time direction */
	if (!sf_getint("nxw", &nxw) ) nxw=5;
	/* half window size for x2 */
	if (!sf_getbool("twod",&twod) ) twod=true;
	/* y: only twod coherence */
	if (!sf_getbool("verb",&verb) ) verb=true;
	/* verbosity */
	if(sf_getstring("idip")!=NULL) idip  = sf_output("idip");
	/* inline dip */
	else idip = NULL;
	if (!sf_getint("min1", &min1) ) min1=-2;
	if (!sf_getint("max1", &max1) ) max1=2;
	/* inline slope */
	if (!sf_getint("min2", &min2) ) min2=-2;
	if (!sf_getint("max2", &max2) ) max2=2;
	/* xline slope */


	if(n3==1)  twod=true;

	u1 = sf_floatalloc2(n1, n2);
	p1 = sf_floatalloc2(n1, n2);

	h1 = fcoh2_init(n1, n2, min1, max1, ntw);

	if(twod) for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);

		fcoh2_2d(h1, u1, p1, nxw);
		
		sf_floatwrite(u1[0], n1*n2, out);
		if(idip) sf_floatwrite(p1[0], n1*n2, idip);
	}else {  /* end 2D begin 3D */
	
	}
	fcoh2_close(h1);
	free(*u1);
	free(u1);
	free(*p1);
	free(p1);
	return 0;
}






