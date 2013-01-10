/* orientation estimation by structural gradient tensor */

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
#include "rgradient.h"
#include "vecfilt.h"
#include "tls2.h"
#include "runtime.h"

int main(int argc, char*argv[])
{
	sf_file in, wgt, out, az;
	int n1, n2, n3, rect[3], order, n12;
	int i1, i3;
	float **u1, *w, **v1, *v2, a;
	char *interp, *filt;
	void *h1=NULL, *h2=NULL;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (!sf_histint(in, "n3", &n3)) n3=1;

	if (!sf_getint("rect1",&rect[0])) rect[0]=0;
	/* smoothness on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=0;
	/* smoothness on 2nd axis */
	if (!sf_getint("rect3",&rect[2])) rect[2]=0;
	/* smoothness on 3rd axis */
	if (!sf_getint("order",&order)) order=2;
	/* approximating order of finite difference */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */
	if ((filt=sf_getstring("filt"))==NULL) filt="tls";
	/* filter type : ls, tls, tensor */
	if(sf_getstring("weight")!=NULL) wgt  = sf_input("weight");
	else wgt = NULL;
	if(sf_getstring("azimuth")!=NULL) az  = sf_output("azimuth");
	else az = NULL;

	n12 =n1*n2;

	rgradient_init(interp, order, n1, n2);
	if(strcmp(filt,"tensor")==0)
		vecfilt_init(3, n1, n2, rect);
	else if(strcmp(filt,"ls")==0)
	{
		h1 = tls2_init(n1, n2, rect, false, false);
		if(az)
		h2 = tls2_init(n1, n2, rect, false, false);
	}else if(strcmp(filt,"tls")==0)
	{
		h1 = tls2_init(n1, n2, rect, true, false);
		if(az)
		h2 = tls2_init(n1, n2, rect, true, false);
	}else sf_error("filt=%s not support", filt);

	u1 = sf_floatalloc2(n1, n2);
	v1 = sf_floatalloc2(n1*3, n2);
	v2 = sf_floatalloc(n12*3);

	if(wgt) w  = sf_floatalloc(n1*n2);
	else w = NULL;

	runtime_init(n1*n2*sizeof(float));

//	time sequences for recursive operators:
//			0	od	 od+rc		  n3  n3+od	 n3+od+rc
//	i3		|----|-----|-----------|----|-----|
//	read	|----------------------|
//	u1=0						   |----|
//	grad	|----|----------------------|
//	v1=0								|-----|
//	weight		 |----------------------|
//	vecf		 |-----|----------------------|
//	write			   |----------------------|

	for(i3=0; i3<n3+order+rect[2]; i3++)
	{
		if(i3<n3) sf_floatread(u1[0], n1*n2, in);
		if(i3>=n3 && i3<n3+order) memset(u1[0], 0, n1*n2*sizeof(float));
		if(i3<n3+order){
			rgradient(u1, v1);
		}else memset(v1[0], 0, 3*n1*n2*sizeof(float));
		if(i3<order) continue;

		if(wgt && i3<n3+order) sf_floatread(w, n1*n2, wgt);
		if(strcmp(filt,"tensor")==0)
		{
			vecfilt(v1[0], v2, w);
			if(i3<order+rect[2]) continue;

			for(i1=0; i1<n1*n2; i1++)
			{
				a = sqrt(v2[i1*3+1]*v2[i1*3+1]+v2[i1*3+2]*v2[i1*3+2]);
				u1[0][i1] = atan2(a, fabs(v2[i1*3]));
			}
			sf_floatwrite(u1[0], n1*n2, out);
			if(az)
			{
				for(i1=0; i1<n1*n2; i1++)
				{
					if(v2[i1*3+1]<0)
					{
						v2[i1*3+1]=-v2[i1*3+1];
						v2[i1*3+2]=-v2[i1*3+2];
					}
					u1[0][i1] = atan2(v2[i1*3+2], v2[i1*3+1]);
				}
				sf_floatwrite(u1[0], n1*n2, az);
			}
		}else{
			for(i1=0; i1<n1*n2; i1++)
			{
				v2[i1] = sqrt(v1[0][i1*3+1]*v1[0][i1*3+1]
						+ v1[0][i1*3+2]*v1[0][i1*3+2]);
				v2[i1+n12] = fabs(v1[0][i1*3]);
			}
			tls2(h1, v2, v2+n12);
			if(i3 >= order+rect[2])
			sf_floatwrite(v2, n1*n2, out);
			if(az){
				for(i1=0; i1<n1*n2; i1++)
				{
					if(v1[0][i1*3+1] >= 0)
					{
						v2[i1] = v1[0][i1*3+2];
						v2[i1+n12] = v1[0][i1*3+1];
					}else{
						v2[i1] = -v1[0][i1*3+2];
						v2[i1+n12] = -v1[0][i1*3+1];
					}
				}
				tls2(h2, v2, v2+n12);
				if(i3 >= order+rect[2])
				sf_floatwrite(v2, n1*n2, az);
			}
		}
		if(i3>order+rect[2])
		{
			a = runtime(1);
			sf_warning("%d of %d, %f MB/sec;", i3-order-rect[2], n3, a);
		}
	}

	rgradient_close();
	if(strcmp(filt,"tensor")==0)
		vecfilt_close();
	else{
		tls2_close(h1);
		if(az) tls2_close(h2);
	}
	free(u1[0]);
	free(u1);
	free(v1[0]);
	free(v1);
	free(v2);

	if(wgt) free(w);

	return 0;
}



